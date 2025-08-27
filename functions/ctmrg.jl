using ITensors


function ctmrg(T::ITensor, Cₗᵤ::ITensor, Aₗ::ITensor; χmax::Int, cutoff=0.0, nsteps::Int)
    """
    this is ISOTROPIC ctmrg 
    Make sure T has indices sh, sh', sv, sv' 
            Aₗ  has indices (lᵥ, lᵥ', sₕ)
            Cₗᵤ has indices (lᵥ, lₕ)
    E.G.
    T = aklt_norm_tensor([sₕ', sᵥ', sₕ, sᵥ])
    χ0 = 1
    l = Index(χ0, "Link")
    lₕ = addtags(l, "horiz")
    lᵥ = addtags(l, "vert")
    # Initial CTM
    Cₗᵤ = ITensor(lᵥ, lₕ)
    Cₗᵤ[1, 1] = 1.0
    # Initial HRTM
    Aₗ = ITensor(lᵥ, lᵥ', sₕ)
    Aₗ[lᵥ => 1, lᵥ' => 1, sₕ => 1] = 1.0
    Cₗᵤ, Aₗ = ctmrg(T, Cₗᵤ, Aₗ; χmax=χmax, cutoff=cutoff, nsteps=nsteps)
    """
  sₕ = commonind(T, Aₗ)
  sᵥ = uniqueind(T, Aₗ, Aₗ'; plev=0)
  lᵥ = commonind(Cₗᵤ, Aₗ)
  lₕ = uniqueind(Cₗᵤ, Aₗ)
  Aᵤ = replaceinds(Aₗ, lᵥ => lₕ, lᵥ' => lₕ', sₕ => sᵥ)
  Cₗᵤ = dense(Cₗᵤ)
  for i in 1:nsteps
    ## Get the grown corner transfer matrix (CTM)
    Cₗᵤ⁽¹⁾ = Aₗ * Cₗᵤ * Aᵤ * T

    ## Diagonalize the grown CTM
    # TODO: replace with
    # eigen(Cₗᵤ⁽¹⁾, "horiz" => "vert"; tags = "horiz" => "vert", kwargs...)
    Cₗᵤ, Uᵥ = eigen(
      Cₗᵤ⁽¹⁾,
      (lₕ', sₕ'),
      (lᵥ', sᵥ');
      ishermitian=true,
      cutoff,
      maxdim=χmax,
      lefttags=tags(lₕ),
      righttags=tags(lᵥ),
    )
    Cₗᵤ = dense(Cₗᵤ)
    lᵥ = commonind(Cₗᵤ, Uᵥ)
    lₕ = uniqueind(Cₗᵤ, Uᵥ)

    # The renormalized CTM is the diagonal matrix of eigenvalues
    # Normalize the CTM
    Cₗ = Cₗᵤ * prime(dag(Cₗᵤ), lₕ)
    normC = (Cₗ * dag(Cₗ))[]^(1 / 4)
    Cₗᵤ = Cₗᵤ / normC

    # Calculate the renormalized half row transfer matrix (HRTM)
    Uᵥ = noprime(Uᵥ)
    Aₗ = Aₗ * Uᵥ * T * dag(Uᵥ')
    Aₗ = replaceinds(Aₗ, sₕ' => sₕ)

    # Normalize the HRTM
    ACₗ = Aₗ * Cₗᵤ * prime(dag(Cₗᵤ))
    normA = √((ACₗ * dag(ACₗ))[])
    Aₗ = Aₗ / normA
    Aᵤ = replaceinds(Aₗ, lᵥ => lₕ, lᵥ' => lₕ', sₕ => sᵥ)
  end
  return Cₗᵤ, Aₗ
end