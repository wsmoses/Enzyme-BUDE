include("BUDE.jl")
using StaticArrays
using Enzyme
import Base.Experimental: @aliasscope

Enzyme.API.printperf!(true)

const Device = (undef, "CPU", "Threaded")

function devices()
  return [Device]
end

const Float32Max = typemax(Float32)

function run(params::Params, deck::Deck, _::DeviceWithRepr)
  println("Using max $(Threads.nthreads()) threads")
  if params.ppwi != DefaultPPWI
    @warn "Threaded implementation only uses wgsize, the PPWI argument is ignored"
  end

  poses = size(deck.poses)[2]
  etotals = Array{Float32}(undef, poses)

  d_protein = fill!(similar(deck.protein), Atom(0f0, 0f0, 0f0, 0))
  d_etotals = fill!(similar(etotals), 0)
  # warmup
  if params.enzyme
     @time autodiff(fasten_main, Const,
        Val(convert(Int, params.wgsize)),
        Duplicated(deck.protein, d_protein),
        Const(deck.ligand),
        Const(deck.forcefield),
        Const(deck.poses),
        Duplicated(etotals, d_etotals)
      )
  else
    @time fasten_main(
      Val(convert(Int, params.wgsize)),
      deck.protein,
      deck.ligand,
      deck.forcefield,
      deck.poses,
      etotals,
    )
  end

  if params.enzyme
    elapsed = @elapsed for _ = 1:params.iterations
      autodiff(fasten_main, Const,
        Val(convert(Int, params.wgsize)),
        Duplicated(deck.protein, d_protein),
        Const(deck.ligand),
        Const(deck.forcefield),
        Const(deck.poses),
        Duplicated(etotals, d_etotals)
      )
    end
  else
    elapsed = @elapsed for _ = 1:params.iterations
      fasten_main(
        Val(convert(Int, params.wgsize)),
        deck.protein,
        deck.ligand,
        deck.forcefield,
        deck.poses,
        etotals,
      )
    end
  end
  
  (etotals, elapsed, params.wgsize)
end

@eval @fastmath function fasten_main(
  ::Val{WGSIZE},
  protein::AbstractArray{Atom},
  ligand::AbstractArray{Atom},
  forcefield::AbstractArray{FFParams},
  poses::AbstractArray{Float32,2},
  etotals::AbstractArray{Float32}, # Only output
) where {WGSIZE}
  protein = Base.Experimental.Const(protein)
  ligand = Base.Experimental.Const(ligand)
  forcefield = Base.Experimental.Const(forcefield)
  poses = Base.Experimental.Const(poses)
  nposes::Int = size(poses)[2]
  numGroups::Int = nposes ÷ WGSIZE

  function inner(group)
    nligand::Int = length(ligand)
    nprotein::Int = length(protein)
    etot = MArray{Tuple{WGSIZE}, Float32}(undef)
    transform = MArray{Tuple{WGSIZE, 3, 4},Float32}(undef)

    @aliasscope @simd ivdep for i = 1:WGSIZE
    # for i = 1:WGSIZE
      ix = (group - 1) * (WGSIZE) + i
      @inbounds sx::Float32 = sin(poses[1, ix])
      @inbounds cx::Float32 = cos(poses[1, ix])
      @inbounds sy::Float32 = sin(poses[2, ix])
      @inbounds cy::Float32 = cos(poses[2, ix])
      @inbounds sz::Float32 = sin(poses[3, ix])
      @inbounds cz::Float32 = cos(poses[3, ix])
      @inbounds transform[i, 1, 1] = cy * cz
      @inbounds transform[i, 1, 2] = sx * sy * cz - cx * sz
      @inbounds transform[i, 1, 3] = cx * sy * cz + sx * sz
      @inbounds transform[i, 1, 4] = poses[4, ix]
      @inbounds transform[i, 2, 1] = cy * sz
      @inbounds transform[i, 2, 2] = sx * sy * sz + cx * cz
      @inbounds transform[i, 2, 3] = cx * sy * sz - sx * cz
      @inbounds transform[i, 2, 4] = poses[5, ix]
      @inbounds transform[i, 3, 1] = -sy
      @inbounds transform[i, 3, 2] = sx * cy
      @inbounds transform[i, 3, 3] = cx * cy
      @inbounds transform[i, 3, 4] = poses[6, ix]
      @inbounds etot[i] = Zero
      # $(Expr(:loopinfo, Symbol("julia.simdloop"), (Symbol("llvm.loop.unroll.full"),)))
    end

    function inner_ligand(il, transform)
      l_etot = zero(MArray{Tuple{WGSIZE}, Float32})

      @inbounds l_atom::Atom = ligand[il]

      @inbounds l_params::FFParams = forcefield[l_atom.type+1]
      lhphb_ltz::Bool = l_params.hphb < Zero
      lhphb_gtz::Bool = l_params.hphb > Zero

      lpos = MArray{Tuple{WGSIZE, 3}, Float32}(undef)

      @simd ivdep for i = 1:WGSIZE
      # for i = 1:WGSIZE
        @inbounds lpos[i, 1] = (
          transform[i, 1, 4] +
          l_atom.x * transform[i, 1, 1] +
          l_atom.y * transform[i, 1, 2] +
          l_atom.z * transform[i, 1, 3]
        )
        @inbounds lpos[i, 2] = (
          transform[i, 2, 4] +
          l_atom.x * transform[i, 2, 1] +
          l_atom.y * transform[i, 2, 2] +
          l_atom.z * transform[i, 2, 3]
        )
        @inbounds lpos[i, 3] = (
          transform[i, 3, 4] +
          l_atom.x * transform[i, 3, 1] +
          l_atom.y * transform[i, 3, 2] +
          l_atom.z * transform[i, 3, 3]
        )
        # $(Expr(:loopinfo, Symbol("julia.simdloop"), (Symbol("llvm.loop.unroll.full"),)))
      end

      @aliasscope for ip = 1:nprotein
        @inbounds p_atom = protein[ip]
        @inbounds p_params::FFParams = forcefield[p_atom.type+1]

        radij::Float32 = p_params.radius + l_params.radius
        r_radij::Float32 = One / radij

        elcdst::Float32 = (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Four : Two
        elcdst1::Float32 =
          (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Quarter : Half
        type_E::Bool = ((p_params.hbtype == HbtypeE || l_params.hbtype == HbtypeE))

        phphb_ltz::Bool = p_params.hphb < Zero
        phphb_gtz::Bool = p_params.hphb > Zero
        phphb_nz::Bool = p_params.hphb != Zero
        p_hphb::Float32 = p_params.hphb * (phphb_ltz && lhphb_gtz ? -One : One)
        l_hphb::Float32 = l_params.hphb * (phphb_gtz && lhphb_ltz ? -One : One)
        distdslv::Float32 =
          (phphb_ltz ? (lhphb_ltz ? Npnpdist : Nppdist) : (lhphb_ltz ? Nppdist : -Float32Max))
        r_distdslv::Float32 = One / distdslv

        chrg_init::Float32 = l_params.elsc * p_params.elsc
        dslv_init::Float32 = p_hphb + l_hphb

        @simd ivdep for i = 1:WGSIZE
        # for i = 1:WGSIZE
          @inbounds x::Float32 = lpos[i, 1] - p_atom.x
          @inbounds y::Float32 = lpos[i, 2] - p_atom.y
          @inbounds z::Float32 = lpos[i, 3] - p_atom.z

          distij::Float32 = sqrt(x ^ 2 + y ^ 2+ z ^ 2)

          # Calculate the sum of the sphere radii
          distbb::Float32 = distij - radij
          zone1::Bool = (distbb < Zero)

          # Calculate steric energy
          @inbounds l_etot[i] += (One - (distij * r_radij)) * (zone1 ? Two * Hardness : Zero)

          # Calculate formal and dipole charge interactions
          chrg_e::Float32 =
            chrg_init * ((zone1 ? One : (One - distbb * elcdst1)) * (distbb < elcdst ? One : Zero))
          neg_chrg_e::Float32 = -abs(chrg_e)
          chrg_e = type_E ? neg_chrg_e : chrg_e
          @inbounds l_etot[i] += chrg_e * Cnstnt

          # Calculate the two cases for Nonpolar-Polar repulsive interactions
          coeff::Float32 = (One - (distbb * r_distdslv))
          dslv_e::Float32 = dslv_init * ((distbb < distdslv && phphb_nz) ? One : Zero)
          dslv_e *= (zone1 ? One : coeff)
          @inbounds l_etot[i] += dslv_e
          # $(Expr(:loopinfo, Symbol("julia.simdloop"), (Symbol("llvm.loop.unroll.full"),)))
        end
      end
      return SArray(l_etot)
    end
    for il = 1:nligand
      etot .+= inner_ligand(il, SArray(transform))
    end
    @simd ivdep for i = 1:WGSIZE
    # for i = 1:WGSIZE
      ix = (group - 1) * WGSIZE + i
      @inbounds etotals[ix] = etot[i] * Half
      # $(Expr(:loopinfo, Symbol("julia.simdloop"), (Symbol("llvm.loop.unroll.full"),)))
    end
  end

  Threads.@threads for group = 1:numGroups
    inner(group)
  end
end

main()
