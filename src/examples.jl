include("Polyhedron.jl")
include("SimplicialSurface.jl")
include("merging.jl")
include("Combinatorics//combinatorics.jl")
include("affine_geometry.jl")

##################################################################################
################# waterbomb
##################################################################################
function waterbomb_cells(eta::Real, zeta::Real, gamma::Real, del_beta::Real=0, base_length::Real=1)
  # translate the following code from Python to Julia
  a = base_length
  b = eta * a
  c = zeta * a
  beta_sym = acos((a * (1 - 2 * sin(gamma))) / sqrt(a^2 + b^2))
  beta = beta_sym + del_beta


  # Koordinaten siehe Maple Sheet
  O = [0, 0, 0]
  LOparam_simpl = [(-a * sin(gamma) * sin(2 * gamma) - cos(gamma) * cos(2 * gamma) * a + cos(gamma) * sqrt(a^2 + b^2) * cos(beta)) / sin(2 * gamma),
    sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) / sin(2 * gamma),
    (a * cos(gamma) * sin(2 * gamma) - sin(gamma) * cos(2 * gamma) * a + sin(gamma) * sqrt(a^2 + b^2) * cos(beta)) / sin(2 * gamma)]
  ROparam_simpl = [(sqrt((-2 * a^3 + 2 * a * b^2) * cos(beta) * sqrt(a^2 + b^2) - (a^2 + b^2)^2 * cos(beta)^2 - a^4 + 3 * a^2 * b^2) * sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) * cos(gamma) - cos(gamma) * cos(beta) * ((a^2 - b^2) * cos(2 * gamma) - a^2) * sqrt(a^2 + b^2) - (sin(gamma) * (cos(beta) - 1) * (cos(beta) + 1) * (a^2 + b^2) * sin(2 * gamma) + (cos(beta)^2 * (a^2 + b^2) * cos(2 * gamma) - a^2 + b^2) * cos(gamma)) * a) / (a^2 + b^2) / sin(beta)^2 / sin(2 * gamma),
    ((b^2 - a^2 - a * sqrt(a^2 + b^2) * cos(beta)) * sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) + sqrt((-2 * a^3 + 2 * a * b^2) * cos(beta) * sqrt(a^2 + b^2) - (a^2 + b^2)^2 * cos(beta)^2 - a^4 + 3 * a^2 * b^2) * (a - sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma))) / (a^2 + b^2) / sin(beta)^2 / sin(2 * gamma),
    (-sqrt((-2 * a^3 + 2 * a * b^2) * cos(beta) * sqrt(a^2 + b^2) - (a^2 + b^2)^2 * cos(beta)^2 - a^4 + 3 * a^2 * b^2) * sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) * sin(gamma) + cos(beta) * ((a^2 - b^2) * cos(2 * gamma) - a^2) * sin(gamma) * sqrt(a^2 + b^2) - (cos(gamma) * (cos(beta) - 1) * (cos(beta) + 1) * (a^2 + b^2) * sin(2 * gamma) - (cos(beta)^2 * (a^2 + b^2) * cos(2 * gamma) - a^2 + b^2) * sin(gamma)) * a) / (a^2 + b^2) / sin(beta)^2 / sin(2 * gamma)]
  MR_simpl = [c * sin(gamma), 0, c * cos(gamma)]
  RUparam_simpl = [(a * sin(gamma) * sin(2 * gamma) + cos(gamma) * cos(2 * gamma) * a - cos(gamma) * sqrt(a^2 + b^2) * cos(beta)) / sin(2 * gamma),
    -sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) / sin(2 * gamma),
    (a * cos(gamma) * sin(2 * gamma) - sin(gamma) * cos(2 * gamma) * a + sin(gamma) * sqrt(a^2 + b^2) * cos(beta)) / sin(2 * gamma)]
  LUparam_simpl = [(-sqrt((-2 * a^3 + 2 * a * b^2) * cos(beta) * sqrt(a^2 + b^2) - (a^2 + b^2)^2 * cos(beta)^2 - a^4 + 3 * a^2 * b^2) * sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) * cos(gamma) + cos(gamma) * cos(beta) * ((a^2 - b^2) * cos(2 * gamma) - a^2) * sqrt(a^2 + b^2) + (sin(gamma) * (cos(beta) - 1) * (cos(beta) + 1) * (a^2 + b^2) * sin(2 * gamma) + (cos(beta)^2 * (a^2 + b^2) * cos(2 * gamma) - a^2 + b^2) * cos(gamma)) * a) / (a^2 + b^2) / sin(beta)^2 / sin(2 * gamma),
    ((a * sqrt(a^2 + b^2) * cos(beta) + a^2 - b^2) * sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) - sqrt((-2 * a^3 + 2 * a * b^2) * cos(beta) * sqrt(a^2 + b^2) - (a^2 + b^2)^2 * cos(beta)^2 - a^4 + 3 * a^2 * b^2) * (a - sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma))) / (a^2 + b^2) / sin(beta)^2 / sin(2 * gamma),
    (-sqrt((-2 * a^3 + 2 * a * b^2) * cos(beta) * sqrt(a^2 + b^2) - (a^2 + b^2)^2 * cos(beta)^2 - a^4 + 3 * a^2 * b^2) * sqrt(2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) + (-a^2 - b^2) * cos(2 * gamma)^2 + (-a^2 - b^2) * cos(beta)^2 + b^2) * sin(gamma) + cos(beta) * ((a^2 - b^2) * cos(2 * gamma) - a^2) * sin(gamma) * sqrt(a^2 + b^2) - (cos(gamma) * (cos(beta) - 1) * (cos(beta) + 1) * (a^2 + b^2) * sin(2 * gamma) - (cos(beta)^2 * (a^2 + b^2) * cos(2 * gamma) - a^2 + b^2) * sin(gamma)) * a) / (a^2 + b^2) / sin(beta)^2 / sin(2 * gamma)]
  ML_simpl = [-c * sin(gamma), 0, c * cos(gamma)]

  cos_psi1 = ((b^2 - a^2) - a * sqrt(a^2 + b^2) * cos(beta)) / (b * sqrt(a^2 + b^2) * sin(beta))
  sin_psi1 = sqrt(a^2 * (3 * b^2 - a^2) + 2 * a * (b^2 - a^2) * sqrt(a^2 + b^2) * cos(beta) - (a^2 + b^2)^2 * cos(beta)^2) / (b * sqrt(a^2 + b^2) * sin(beta))
  cos_psi5 = (sqrt(a^2 + b^2) * cos(beta) - a * cos(2 * gamma)) / (b * sin(2 * gamma))
  sin_psi5 = sqrt(b^2 + 2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) - (a^2 + b^2) * (cos(beta)^2 + cos(2 * gamma)^2)) / (b * sin(2 * gamma))
  cos_psi6 = (a - sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma)) / (sqrt(a^2 + b^2) * sin(beta) * sin(2 * gamma))
  sin_psi6 = sqrt(b^2 + 2 * a * sqrt(a^2 + b^2) * cos(beta) * cos(2 * gamma) - (a^2 + b^2) * (cos(beta)^2 + cos(2 * gamma)^2)) / (sqrt(a^2 + b^2) * sin(beta) * sin(2 * gamma))

  cos_phi1 = cos_psi1 * cos_psi6 - sin_psi1 * sin_psi6
  sin_phi1 = sin_psi1 * cos_psi6 + cos_psi1 * sin_psi6
  cos_phi2 = cos_psi5
  sin_phi2 = sin_psi5

  sin_phi1phi2 = sin_psi1 * cos_psi6 * cos_psi5 - sin_psi1 * sin_psi6 * sin_psi5 + cos_psi1 * sin_psi6 * cos_psi5 + cos_psi1 * cos_psi6 * sin_psi5
  cos_phi1phi2 = cos_psi1 * cos_psi6 * cos_psi5 - cos_psi1 * sin_psi6 * sin_psi5 - sin_psi1 * sin_psi5 * cos_psi5 - sin_psi1 * cos_psi6 * sin_psi5

  W = 2 * a * (1 - cos_phi1phi2) * (sin(2 * gamma)^2 * a * c^2 - (cos_phi2 + cos_phi1) * (cos(2 * gamma) - 1) * sin(2 * gamma) * b * c^2 - 2 * (cos_phi2 + cos_phi1) * sin(2 * gamma) * a * b * c - 2 * a * c * (2 * a - c) * (cos(2 * gamma) + 1) + 2 * a * b^2 * (cos_phi1phi2 + 1)) - (((a - c)^2 + b^2) * (cos_phi2^2 + cos_phi1^2) - 2 * cos_phi1 * cos_phi2 * ((a - c)^2 + cos_phi1phi2 * b^2)) * c^2 * sin(2 * gamma)^2
  T_P = ((a - c)^2 + b^2) * (cos_phi1phi2 - 1) * (a * (2 * b^2 * (cos_phi1phi2 + 1) + 4 * a * (a - c)) * cos(2 * gamma) + 2 * b * sin(2 * gamma) * (b^2 * cos_phi2 * (cos_phi1phi2 + 1) + (a - c) * ((a - c) * cos_phi2 + a * cos_phi2 + c * cos_phi1)) + 2 * a * (2 * a * (a - c) - b^2 * (cos_phi1phi2 + 1)))
  # t = (2*b*(a-c)*(cos_phi1phi2 - 1)*sqrt((a-c)^2 + b^2)*sqrt(W) - 2*b*(a*b*c*sin_phi1phi2*(cos_phi1phi2 - 1)*cos(2*gamma) + (((a-c)^2 + b^2*cos_phi1phi2)*cos_phi2 - cos_phi1*((a-c)^2 + b^2))*sin_phi1phi2*c*sin(2*gamma) + a*b*sin_phi1phi2*(cos_phi1phi2 - 1)*(2*a - c))*sqrt((a-c)^2 + b^2)+ (cos_phi1phi2 - 1)*(2*b^2*cos_phi1phi2 + 4*(a-c)^2 + 2*b^2)*sqrt((a-c)^2 + b^2)*sqrt(4*a*b^2*(cos_phi2*sin(2*gamma)*b + (cos(2*gamma) + 1)*a) - (a^2 + b^2)*(cos_phi2*sin(2*gamma)*b + (cos(2*gamma) + 1)*a)^2)) / ((2*b*((a-c)^2 + b^2)*sin_phi1phi2 * sqrt(W) - T_P))
  t = (2 * b * (a - c) * (cos_phi1phi2 - 1) * sqrt((a - c)^2 + b^2) * sqrt(W) - 2 * b * (a * b * c * sin_phi1phi2 * (cos_phi1phi2 - 1) * cos(2 * gamma) + (((a - c)^2 + b^2 * cos_phi1phi2) * cos_phi2 - cos_phi1 * ((a - c)^2 + b^2)) * sin_phi1phi2 * c * sin(2 * gamma) + a * b * sin_phi1phi2 * (cos_phi1phi2 - 1) * (2 * a - c)) * sqrt((a - c)^2 + b^2) - (cos_phi1phi2 - 1) * (2 * b^2 * cos_phi1phi2 + 4 * (a - c)^2 + 2 * b^2) * sqrt((a - c)^2 + b^2) * sqrt(4 * a * b^2 * (cos_phi2 * sin(2 * gamma) * b + (cos(2 * gamma) + 1) * a) - (a^2 + b^2) * (cos_phi2 * sin(2 * gamma) * b + (cos(2 * gamma) + 1) * a)^2)) / ((2 * b * ((a - c)^2 + b^2) * sin_phi1phi2 * sqrt(W) - T_P))
  # t = (    -2*b*(a-c)*(cos_phi1phi2 - 1)*sqrt((a-c)^2 + b^2)*sqrt(W) - 2*b*(a*b*c*sin_phi1phi2*(cos_phi1phi2 - 1)*cos(2*gamma) + (((a-c)^2 + b^2*cos_phi1phi2)*cos_phi2 - cos_phi1*((a-c)^2 + b^2))*sin_phi1phi2*c*sin(2*gamma) + a*b*sin_phi1phi2*(cos_phi1phi2 - 1)*(2*a - c))*sqrt((a-c)^2 + b^2)+ (cos_phi1phi2 - 1)*(2*b^2*cos_phi1phi2 + 4*(a-c)^2 + 2*b^2)*sqrt((a-c)^2 + b^2)*sqrt(4*a*b^2*(cos_phi2*sin(2*gamma)*b + (cos(2*gamma) + 1)*a) - (a^2 + b^2)*(cos_phi2*sin(2*gamma)*b + (cos(2*gamma) + 1)*a)^2)) / ((-2*b*((a-c)^2 + b^2)*sin_phi1phi2 * sqrt(W) - T_P))
  # t = (    -2*b*(a-c)*(cos_phi1phi2 - 1)*sqrt((a-c)^2 + b^2)*sqrt(W) - 2*b*(a*b*c*sin_phi1phi2*(cos_phi1phi2 - 1)*cos(2*gamma) + (((a-c)^2 + b^2*cos_phi1phi2)*cos_phi2 - cos_phi1*((a-c)^2 + b^2))*sin_phi1phi2*c*sin(2*gamma) + a*b*sin_phi1phi2*(cos_phi1phi2 - 1)*(2*a - c))*sqrt((a-c)^2 + b^2)- (cos_phi1phi2 - 1)*(2*b^2*cos_phi1phi2 + 4*(a-c)^2 + 2*b^2)*sqrt((a-c)^2 + b^2)*sqrt(4*a*b^2*(cos_phi2*sin(2*gamma)*b + (cos(2*gamma) + 1)*a) - (a^2 + b^2)*(cos_phi2*sin(2*gamma)*b + (cos(2*gamma) + 1)*a)^2)) / ((-2*b*((a-c)^2 + b^2)*sin_phi1phi2 * sqrt(W) - T_P))


  T_Q = ((a - c)^2 + b^2) * (cos_phi1phi2 - 1) * (a * (2 * b^2 * (cos_phi1phi2 + 1) + 4 * a * (a - c)) * cos(2 * gamma) + 2 * b * sin(2 * gamma) * (b^2 * cos_phi1 * (cos_phi1phi2 + 1) + (a - c) * ((a - c) * cos_phi1 + a * cos_phi1 + c * cos_phi2)) + 2 * a * (2 * a * (a - c) - b^2 * (cos_phi1phi2 + 1)))
  # u = (2*b*(a-c)*(cos_phi1phi2 - 1)*sqrt((a-c)^2 + b^2)*sqrt(W) - 2*b*(a*b*c*sin_phi1phi2*(cos_phi1phi2 - 1)*cos(2*gamma) + (((a-c)^2 + b^2*cos_phi1phi2)*cos_phi1 - cos_phi2*((a-c)^2 + b^2))*sin_phi1phi2*c*sin(2*gamma) + a*b*sin_phi1phi2*(cos_phi1phi2 - 1)*(2*a - c))*sqrt((a-c)^2 + b^2)+ (cos_phi1phi2 - 1)*(2*b^2*cos_phi1phi2 + 4*(a-c)^2 + 2*b^2)*sqrt((a-c)^2 + b^2)*sqrt(4*a*b^2*(cos_phi1*sin(2*gamma)*b + (cos(2*gamma) + 1)*a) - (a^2 + b^2)*(cos_phi1*sin(2*gamma)*b + (cos(2*gamma) + 1)*a)^2)) / ((2*b*((a-c)^2 + b^2)*sin_phi1phi2 * sqrt(W) - T_Q))
  u = (2 * b * (a - c) * (cos_phi1phi2 - 1) * sqrt((a - c)^2 + b^2) * sqrt(W) - 2 * b * (a * b * c * sin_phi1phi2 * (cos_phi1phi2 - 1) * cos(2 * gamma) + (((a - c)^2 + b^2 * cos_phi1phi2) * cos_phi1 - cos_phi2 * ((a - c)^2 + b^2)) * sin_phi1phi2 * c * sin(2 * gamma) + a * b * sin_phi1phi2 * (cos_phi1phi2 - 1) * (2 * a - c)) * sqrt((a - c)^2 + b^2) - (cos_phi1phi2 - 1) * (2 * b^2 * cos_phi1phi2 + 4 * (a - c)^2 + 2 * b^2) * sqrt((a - c)^2 + b^2) * sqrt(4 * a * b^2 * (cos_phi1 * sin(2 * gamma) * b + (cos(2 * gamma) + 1) * a) - (a^2 + b^2) * (cos_phi1 * sin(2 * gamma) * b + (cos(2 * gamma) + 1) * a)^2)) / ((2 * b * ((a - c)^2 + b^2) * sin_phi1phi2 * sqrt(W) - T_Q))
  # u = (    -2*b*(a-c)*(cos_phi1phi2 - 1)*sqrt((a-c)^2 + b^2)*sqrt(W) - 2*b*(a*b*c*sin_phi1phi2*(cos_phi1phi2 - 1)*cos(2*gamma) + (((a-c)^2 + b^2*cos_phi1phi2)*cos_phi1 - cos_phi2*((a-c)^2 + b^2))*sin_phi1phi2*c*sin(2*gamma) + a*b*sin_phi1phi2*(cos_phi1phi2 - 1)*(2*a - c))*sqrt((a-c)^2 + b^2)+ (cos_phi1phi2 - 1)*(2*b^2*cos_phi1phi2 + 4*(a-c)^2 + 2*b^2)*sqrt((a-c)^2 + b^2)*sqrt(4*a*b^2*(cos_phi1*sin(2*gamma)*b + (cos(2*gamma) + 1)*a) - (a^2 + b^2)*(cos_phi1*sin(2*gamma)*b + (cos(2*gamma) + 1)*a)^2)) / ((-2*b*((a-c)^2 + b^2)*sin_phi1phi2 * sqrt(W) - T_Q))
  # u = (    -2*b*(a-c)*(cos_phi1phi2 - 1)*sqrt((a-c)^2 + b^2)*sqrt(W) - 2*b*(a*b*c*sin_phi1phi2*(cos_phi1phi2 - 1)*cos(2*gamma) + (((a-c)^2 + b^2*cos_phi1phi2)*cos_phi1 - cos_phi2*((a-c)^2 + b^2))*sin_phi1phi2*c*sin(2*gamma) + a*b*sin_phi1phi2*(cos_phi1phi2 - 1)*(2*a - c))*sqrt((a-c)^2 + b^2)- (cos_phi1phi2 - 1)*(2*b^2*cos_phi1phi2 + 4*(a-c)^2 + 2*b^2)*sqrt((a-c)^2 + b^2)*sqrt(4*a*b^2*(cos_phi1*sin(2*gamma)*b + (cos(2*gamma) + 1)*a) - (a^2 + b^2)*(cos_phi1*sin(2*gamma)*b + (cos(2*gamma) + 1)*a)^2)) / ((-2*b*((a-c)^2 + b^2)*sin_phi1phi2 * sqrt(W) - T_Q))

  P = [(2 * sqrt(a^2 - 2 * a * c + b^2 + c^2) * cos(gamma) * sin_phi1 * b * c * t + ((t^2 + 1) * c^3 + (-t^2 - 1) * a * c^2 + ((t^2 - 1) * b^2 - a^2 * (t^2 + 1)) * c + a * (t^2 + 1) * (a^2 + b^2)) * sin(gamma) - ((t^2 - 1) * c^2 - 2 * a * c * t^2 + (t^2 + 1) * (a^2 + b^2)) * cos(gamma) * cos_phi1 * b) / (t^2 + 1) / (a^2 - 2 * a * c + b^2 + c^2), (2 * sqrt(a^2 - 2 * a * c + b^2 + c^2) * cos_phi1 * c * t + ((t^2 - 1) * c^2 - 2 * a * c * t^2 + (t^2 + 1) * (a^2 + b^2)) * sin_phi1) * b / (t^2 + 1) / (a^2 - 2 * a * c + b^2 + c^2), (-2 * sqrt(a^2 - 2 * a * c + b^2 + c^2) * sin(gamma) * sin_phi1 * b * c * t + ((t^2 + 1) * c^3 + (-t^2 - 1) * a * c^2 + ((t^2 - 1) * b^2 - a^2 * (t^2 + 1)) * c + a * (t^2 + 1) * (a^2 + b^2)) * cos(gamma) + ((t^2 - 1) * c^2 - 2 * a * c * t^2 + (t^2 + 1) * (a^2 + b^2)) * sin(gamma) * cos_phi1 * b) / (t^2 + 1) / (a^2 - 2 * a * c + b^2 + c^2)]
  Q = [(2 * cos(gamma) * sqrt(a^2 - 2 * a * c + b^2 + c^2) * sin_phi2 * b * c * u + ((u^2 + 1) * c^3 + (-u^2 - 1) * a * c^2 + ((u^2 - 1) * b^2 - a^2 * (u^2 + 1)) * c + a * (u^2 + 1) * (a^2 + b^2)) * sin(gamma) - ((u^2 - 1) * c^2 - 2 * a * c * u^2 + (u^2 + 1) * (a^2 + b^2)) * cos(gamma) * cos_phi2 * b) / (u^2 + 1) / (a^2 - 2 * a * c + b^2 + c^2), -(2 * sqrt(a^2 - 2 * a * c + b^2 + c^2) * cos_phi2 * c * u + ((u^2 - 1) * c^2 - 2 * a * c * u^2 + (u^2 + 1) * (a^2 + b^2)) * sin_phi2) * b / (u^2 + 1) / (a^2 - 2 * a * c + b^2 + c^2), (-2 * sin(gamma) * sqrt(a^2 - 2 * a * c + b^2 + c^2) * sin_phi2 * b * c * u + ((u^2 + 1) * c^3 + (-u^2 - 1) * a * c^2 + ((u^2 - 1) * b^2 - a^2 * (u^2 + 1)) * c + a * (u^2 + 1) * (a^2 + b^2)) * cos(gamma) + ((u^2 - 1) * c^2 - 2 * a * c * u^2 + (u^2 + 1) * (a^2 + b^2)) * sin(gamma) * cos_phi2 * b) / (u^2 + 1) / (a^2 - 2 * a * c + b^2 + c^2)]


  cell1 = hcat(O, LOparam_simpl, ROparam_simpl, MR_simpl, RUparam_simpl, LUparam_simpl, ML_simpl)

  # second cell
  A = affinemap(hcat(O, ML_simpl, LUparam_simpl, cross(ML_simpl, LUparam_simpl)), hcat(P, ROparam_simpl, MR_simpl, P + cross(ROparam_simpl - P, MR_simpl - P)))
  cell2 = hcat([A(cell1[:, i]) for i in 1:size(cell1, 2)]...)

  # third cell
  A = affinemap(hcat(O, ML_simpl, LOparam_simpl, cross(ML_simpl, LOparam_simpl)), hcat(Q, RUparam_simpl, MR_simpl, Q + cross(RUparam_simpl - Q, MR_simpl - Q)))
  cell3 = hcat([A(cell1[:, i]) for i in 1:size(cell1, 2)]...)

  faces = [[1, 2, 3], [1, 3, 4], [1, 4, 5], [1, 5, 6], [1, 6, 7], [1, 7, 2], [8, 3, 9], [8, 9, 10], [8, 10, 11], [8, 11, 12], [8, 12, 4], [8, 4, 3], [13, 4, 12], [13, 12, 14], [13, 14, 15], [13, 15, 16], [13, 16, 5], [13, 5, 4]]
  coordinates = hcat([
    cell1[:, 1], cell1[:, 2], cell1[:, 3], cell1[:, 4], cell1[:, 5], cell1[:, 6], cell1[:, 7],
    cell2[:, 1], cell2[:, 2], cell2[:, 3], cell2[:, 4], cell2[:, 5],
    cell3[:, 1], cell3[:, 4], cell3[:, 5], cell3[:, 6]
  ]...)

  # faces = [[0,1,2], [0,2,3], [0,3,4], [0,4,5], [0,5,6], [0,6,1],
  # [7,8,9], [7,9,10], [7,10,11], [7,11,12], [7,12,13], [7,13,8],
  # [14,15,16], [14,16,17], [14,17,18], [14,18,19], [14,19,20], [14,20,15]]
  # coordinates = vstack((cell1, cell2, cell3))

  pattern1 = SimplicialSurface(verts=coordinates, facets=faces)
  # surf = deepcopy(pattern1)
  # surf = merge(pattern1, surf, [11,12,14], [1,6,5], atol=1e-1)
  # surf1 = glue_surfs_along_path(surf, pattern1, [8,9,10,17,21], [5,4,15,14,13], atol=1)
  # surf2 = copy.deepcopy(surf1)
  # surf2 = glue_surfs_along_path(surf1, surf2, [37,36,35,31,30], [28,19,13,14,15], atol=1)
end


# n-gon
function ngon(n::Integer)
  verts = collect(1:n+1)
  edges = vcat([[1, k] for k in 2:n+1], [[mod1(k, n) + 1, mod1(k + 1, n) + 1] for k in 1:n])
  facets = [[1, mod1(k, n) + 1, mod1(k + 1, n) + 1] for k in 1:n]
  return CombSimplicialSurface(verts, edges, facets)
end

function cube_assembly()
  cube1 = Cube

  shift_x = rigidmap([0 1 0 0; 0 0 1 0; 0 0 0 1], [1 2 1 1; 0 0 1 0; 0 0 0 1])
  shift_y = rigidmap([0 1 0 0; 0 0 1 0; 0 0 0 1], [0 1 0 0; 1 1 2 1; 0 0 0 1])

  cube2 = deepcopy(cube1)
  set_verts!(cube2, shift_x(get_verts(cube2)))
  cube3 = deepcopy(cube2)
  set_verts!(cube3, shift_x(get_verts(cube3)))
  cube4 = deepcopy(cube1)
  set_verts!(cube4, shift_y(get_verts(cube4)))
  cube5 = deepcopy(cube2)
  set_verts!(cube5, shift_y(get_verts(cube5)))
  cube6 = deepcopy(cube3)
  set_verts!(cube6, shift_y(get_verts(cube6)))
  cube7 = deepcopy(cube4)
  set_verts!(cube7, shift_y(get_verts(cube7)))
  cube8 = deepcopy(cube5)
  set_verts!(cube8, shift_y(get_verts(cube8)))
  cube9 = deepcopy(cube6)
  set_verts!(cube9, shift_y(get_verts(cube9)))

  return [cube1, cube2, cube3, cube4, cube5, cube6, cube7, cube8, cube9]
end

##################################################################################
################# n-prisms
##################################################################################
"""
    nprism(n::Integer)

TBW
"""
function nprism(n::Integer)
  alpha = 1 / dist([1, 0], [cos(2 * pi / n), sin(2 * pi / n)])
  verts = vcat([[alpha * cos(2 * pi * k / n), alpha * sin(2 * pi * k / n), 0] for k in 0:n-1], [[alpha * cos(2 * pi * k / n), alpha * sin(2 * pi * k / n), 1] for k in 0:n-1])
  facets = vcat([[1:n...]], [[(n+1):(2*n)...]], [[k, k + 1, k + n + 1, k + n] for k in 1:(n-1)], [[n, 1, n + 1, 2 * n]])
  edges = vcat([[k, k + 1] for k in 1:(n-1)], [[n, 1]], [[n + k, n + k + 1] for k in 1:(n-1)], [[2 * n, n + 1]], [[k, k + n] for k in 1:n])
  return Polyhedron(verts, edges, facets)
end

function tiblock(n1::Integer, n2::Integer, nmerges::Integer; atol::Real=1e-8)
  @assert mod(n1, nmerges) == 0 "Number of blocks on the side needs to divide n1."

  if n2 == 3
    sideelem = nprism(3)
    merge!(sideelem, nprism(3), [[3, 1, 6, 4]], [[1, 2, 4, 5]], atol=atol)
    merge!(sideelem, nprism(3), [[2, 3, 5, 6]], [[1, 2, 4, 5]], atol=atol)
  elseif n2 == 4
    sideelem = nprism(4)
    merge!(sideelem, nprism(4), [[2, 3, 6, 7]], [[1, 2, 5, 6]], atol=atol)
    merge!(sideelem, nprism(4), [[4, 1, 8, 5]], [[1, 2, 5, 6]], atol=atol)
  else
    sideelem = nprism(n2)
  end

  block = nprism(n1)
  for k in 3:Int(n1 / nmerges):length(get_facets(nprism(n1)))
    merge!(block, sideelem, [reverse(get_facets(nprism(n1))[k])], [[n2 + 1, 1, 2, n2 + 2]], atol=atol)
  end

  return block
end


##################################################################################
################# paper assemblies
##################################################################################
function assembly1(n1::Int=6, n2::Int=3, nmerges::Int=2, nblocks::Int=7)
  @assert n1 > 3 "not possible for triangle prisma"
  assembly = Polyhedron[]

  block = flattenfacets(tiblock(n1, n2, nmerges))
  theta = -pi / n1
  rotmat = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
  set_verts!(block, rotmat * get_verts(block))
  push!(assembly, block)

  for i in 2:nblocks
    lastblock = assembly[end]
    newblock = deepcopy(lastblock)
    lastcoords = get_verts(lastblock)
    newcoords = get_verts(newblock)

    preim = [newcoords[:, 1], newcoords[:, 2], newcoords[:, 3], newcoords[:, 1] + cross(newcoords[:, 2] - newcoords[:, 1], newcoords[:, 3] - newcoords[:, 1])]
    if mod(i, 2) == 0
      im = [lastcoords[:, n1+2], lastcoords[:, n1+3], lastcoords[:, n1+4], lastcoords[:, n1+2] + cross(lastcoords[:, n1+3] - lastcoords[:, n1+2], lastcoords[:, n1+4] - lastcoords[:, n1+2])]
    else
      im = [lastcoords[:, 2*n1], lastcoords[:, n1+1], lastcoords[:, n1+2], lastcoords[:, 2*n1] + cross(lastcoords[:, n1+1] - lastcoords[:, 2*n1], lastcoords[:, n1+2] - lastcoords[:, 2*n1])]
    end
    aff = rigidmap(preim, im)
    set_verts!(newblock, aff(newcoords))

    push!(assembly, newblock)
  end

  frame = [1, nblocks]

  return assembly, frame
end

function assembly2(n1::Int=6, n2::Int=3, nmerges::Int=2, nblocks::Int=7)
  @assert n1 > 3 "not possible for triangle prisma"
  assembly = Polyhedron[]

  block = flattenfacets(tiblock(n1, n2, nmerges))
  theta = -pi / n1
  rotmat = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
  set_verts!(block, rotmat * get_verts(block))
  push!(assembly, block)

  for i in 2:nblocks
    lastblock = assembly[end]
    newblock = deepcopy(lastblock)
    lastcoords = get_verts(lastblock)
    newcoords = get_verts(newblock)

    preim = [newcoords[:, 1], newcoords[:, 2], newcoords[:, 3], newcoords[:, 1] + cross(newcoords[:, 2] - newcoords[:, 1], newcoords[:, 3] - newcoords[:, 1])]
    if mod(i, 3) != 0
      im = [lastcoords[:, n1+2], lastcoords[:, n1+3], lastcoords[:, n1+4], lastcoords[:, n1+2] + cross(lastcoords[:, n1+3] - lastcoords[:, n1+2], lastcoords[:, n1+4] - lastcoords[:, n1+2])]
    else
      im = [lastcoords[:, 2*n1-1], lastcoords[:, 2*n1], lastcoords[:, n1+1], lastcoords[:, 2*n1-1] + cross(lastcoords[:, 2*n1] - lastcoords[:, 2*n1-1], lastcoords[:, n1+1] - lastcoords[:, 2*n1-1])]
    end
    aff = rigidmap(preim, im)
    set_verts!(newblock, aff(newcoords))

    push!(assembly, newblock)
  end

  frame = Int[1, nblocks]

  return assembly, frame
end

function assembly3(n1::Int=6, towerheight::Int=7, ntowers::Int=3, nlinks::Int=2; atol::Real=1e-8)
  assembly, frame = assembly2(n1, 3, 2, towerheight)
  link = Int[]

  # offset between blocks of a tower to attach a link between towers
  offset = max(towerheight ÷ 3 ÷ nlinks, 1) * 3
  # corrected number of links that are possible for the towers
  correct_nlinks = length(filter(i -> 4 + (i - 1) * offset <= towerheight, collect(1:nlinks)))
  # base building block
  linkblock = flattenfacets(tiblock(n1, 3, 2))

  function attachlinks(tower)
    for i in 1:correct_nlinks
      newblock = deepcopy(linkblock)
      newcoords = get_verts(newblock)
      towerblocks = tower[[j + (i - 1) * offset for j in collect(1:4)]]
      towerblockcoords = get_verts.(towerblocks)

      preim = newcoords[:, [n1 ÷ 2 + 1, n1 ÷ 2 + 2, n1 ÷ 2 + 1 + n1, n1 ÷ 2 + 2 + n1, collect((2*n1+7):(2*n1+12))...]]
      im = hcat(towerblockcoords[1][:, [2 * n1 + 5, 2 * n1 + 6]], towerblockcoords[4][:, [2 * n1 + 3, 2 * n1 + 4]], towerblockcoords[2][:, [2 * n1, n1 + 1, n1, 1]], towerblockcoords[3][:, [n1 + 2, n1 + 3]])
      display(preim)
      display(im)
      aff = rigidmap(preim, im, atol=atol)
      set_verts!(newblock, aff(newcoords))

      push!(assembly, newblock)
      push!(link, length(assembly))
    end
  end

  attachlinks(assembly)


  for i in 2:ntowers
    newassembly, newframe = assembly2(n1, 3, 2, towerheight)
    newframe = newframe .+ length(assembly)

    append!(frame, newframe)

    # reflink has to attach at the first 4 blocks of the tower
    reflink = assembly[link[end-correct_nlinks+1]]
    refcoords = get_verts(reflink)
    towercoords = get_verts.(newassembly)

    preim = hcat(towercoords[1][:, [2 * n1 + 11, 2 * n1 + 12]], towercoords[4][:, [2 * n1 + 9, 2 * n1 + 10]], towercoords[2][:, [n1 ÷ 2 + n1, n1 ÷ 2 + n1 + 1, n1 ÷ 2, n1 ÷ 2 + 1]], towercoords[3][:, [n1 ÷ 2 + n1 + 2, n1 ÷ 2 + n1 + 3]])
    im = refcoords[:, [1, 2, n1 + 1, n1 + 2, collect((2*n1+1):2*n1+6)]]
    aff = rigidmap(preim, im)

    for block in newassembly
      set_verts!(block, aff(get_verts(block)))
    end
    append!(assembly, newassembly)

    attachlinks(newassembly)
  end

  return assembly, frame, link
end

function assembly4(n1::Int=6, n2::Int=3, nmerges::Int=2, nblocks::Int=7)
  @assert n1 > 3 "not possible for triangle prisma"
  assembly = Polyhedron[]

  block = flattenfacets(tiblock(n1, n2, nmerges))
  theta = -pi / n1
  rotmat = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
  set_verts!(block, [rotmat * coords for coords in get_verts(block)])
  push!(assembly, block)

  for i in 2:nblocks
    lastblock = assembly[end]
    newblock = deepcopy(lastblock)
    lastcoords = get_verts(lastblock)
    newcoords = get_verts(newblock)

    preim = [newcoords[1], newcoords[2], newcoords[3], newcoords[1] + cross(newcoords[2] - newcoords[1], newcoords[3] - newcoords[1])]
    im = [lastcoords[n1+2], lastcoords[n1+3], lastcoords[n1+4], lastcoords[n1+2] + cross(lastcoords[n1+3] - lastcoords[n1+2], lastcoords[n1+4] - lastcoords[n1+2])]
    aff = rigidmap(preim, im)
    set_verts!(newblock, aff(newcoords))

    push!(assembly, newblock)
  end

  frame = Int[1, nblocks]

  return assembly, frame
end

# function tetrahedra_interlocking(n::Int)
#   assembly=[]
#   frame=[]
#   mat=Matrix([[cos(pi/2), sin(pi/2),0.0] [-sin(pi/2),cos(pi/2),0.0] [0.0,0.0,1.0]])
#   vertices=[[0.5,0,0],[-0.5,0.,0.],[0.,0.5,sqrt(3.)/2.],[0.,-0.5,sqrt(3.)/2.]]  
#   for i in 1:n
#     for j in 1:n
#       if i%2==j%2 
#         push!(assembly,map(coor->[(i-1)/2.0,(j-1)/2.0,0.0]+coor,vertices))
#       else
#         push!(assembly,map(coor->[(i-1)/2.0,(j-1)/2.0,0.0]+mat*coor,vertices))
#       end
#       if i==1 || i==n || j==1 || j==n
#         push!(frame,length(assembly))
#       end 
#     end
#   end
#   return map(coor->
#          PolyhedronByVerticesInFacets(coor,[[1,2,3],[1,2,4],[1,3,4],[2,3,4]]),assembly),frame
# end

# function octahedra_interlocking(n::Int)
#   assembly=[]
#   frame=[]
#   vertices=1.0*[[1,0,1],[ 1, 1, 0 ], [ 2, 1, 1 ],[ 1, 1, 2 ], [ 0, 1, 1 ],[ 1, 2, 1 ]]  
#   vec1=vertices[2]-vertices[1]
#   vec2=vertices[3]-vertices[1]
#   for i in 1:n
#     for j in 1:n
#       push!(assembly,map(coor->(i-1)*vec1+(j-1)*vec2+coor,vertices))
#       if i==1 || i==n || j==1 || j==n
#         push!(frame,length(assembly))
#       end 
#     end
#   end
#   temp=[[1,2,3],[1,3,4],[1,4,5],[1,2,5],[6,2,3],[6,3,4],[6,4,5],[6,2,5]]
#   return map(coor->
#          PolyhedronByVerticesInFacets(coor,temp),assembly),frame
# end

# function cube_interlocking(n::Int)
#   assembly=[]
#   frame=[]
#   vertices=1.0*[[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
#   temp=[[1,2,3,4],[5,6,7,8],[5,6,2,1],[6,7,3,2],[7,8,4,3],[1,4,8,5]]
#   vec1=[1.0,-0.5,0.5]
#   vec2=[-0.5,1.0,0.5]
#   for i in 1:n
#     for j in 1:n
#       push!(assembly,map(coor->(i-1)*vec1+(j-1)*vec2+coor,vertices))
#       if i==1 || i==n || j==1 || j==n
#         push!(frame,length(assembly))
#       end 
#     end
#   end
#   return map(coor->
#          PolyhedronByVerticesInFacets(coor,temp),assembly),frame
# end

##################################################################################
################# platonic solids
##################################################################################

####Tetrahedron 

Tetrahedron = Polyhedron([[0.5, 0, 0], [-0.5, 0.0, 0.0], [0.0, 0.5, sqrt(3.0) / 2.0], [0.0, -0.5, sqrt(3.0) / 2.0]],
  [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]],
  [[1, 2, 3], [1, 2, 4], [2, 3, 4], [1, 3, 4]])


####Octahedron
v = 1.0 * [[1, 0, 1], [1, 1, 0], [2, 1, 1], [1, 1, 2], [0, 1, 1], [1, 2, 1]]
f = [[1, 2, 3], [2, 5, 6], [1, 2, 5], [2, 3, 6], [1, 4, 5],
  [3, 4, 6], [1, 3, 4], [4, 5, 6]]
e = [[1, 2], [1, 3], [1, 4], [1, 5], [2, 3], [2, 5], [2, 6],
  [3, 4], [3, 6], [4, 5], [4, 6], [5, 6]]
Octahedron = Polyhedron(v, e, f)

####Cube 
Cube = Polyhedron(1.0 * [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], [[1, 2], [2, 3], [3, 4], [1, 4], [2, 6], [3, 7], [4, 8], [1, 5], [6, 7], [7, 8], [5, 8], [5, 6]], [[1, 2, 3, 4], [5, 6, 7, 8], [5, 6, 2, 1], [6, 7, 3, 2], [7, 8, 4, 3], [1, 4, 8, 5]])

####icosahedron
phi = (1 + sqrt(5.0)) / 2
v = 1.0 * [[0, 1, phi], [0, -1, phi], [0, 1, -phi], [0, -1, -phi], [1, phi, 0], [-1, phi, 0], [1, -phi, 0], [-1, -phi, 0], [phi, 0, 1], [-phi, 0, 1], [phi, 0, -1], [-phi, 0, -1]]

e = [[1, 2], [1, 5], [1, 6], [1, 9], [1, 10], [2, 7], [2, 8],
  [2, 9], [2, 10], [3, 4], [3, 5], [3, 6], [3, 11], [3, 12],
  [4, 7], [4, 8], [4, 11], [4, 12], [5, 6], [5, 9], [5, 11],
  [6, 10], [6, 12], [7, 8], [7, 9], [7, 11], [8, 10], [8, 12],
  [9, 11], [10, 12]]
f = [[1, 2, 9], [1, 2, 10], [1, 5, 6], [1, 5, 9], [1, 6, 10], [2, 7, 8], [2, 7, 9], [2, 8, 10], [3, 4, 11], [3, 4, 12], [3, 5, 6], [3, 5, 11], [3, 6, 12], [4, 7, 8], [4, 7, 11], [4, 8, 12], [5, 9, 11], [6, 10, 12], [7, 9, 11], [8, 10, 12]]
Icosahedron = Polyhedron(v, e, f)

####Dodecahedron
f = [[1, 4, 3, 5, 2], [1, 7, 6, 8, 2], [9, 12, 11, 13, 10], [9, 15, 14, 16, 10], [3, 11, 12, 17, 4], [3, 11, 13, 18, 5], [6, 14, 15, 19, 7], [6, 14, 16, 20, 8], [1, 7, 19, 17, 4], [2, 8, 20, 18, 5], [9, 15, 19, 17, 12], [10, 16, 20, 18, 13]]
v = 1.0 * [[0.539345, 0, 1.41202], [-0.539345, 0, 1.41202], [0, 1.41202, 0.539345], [0.872678, 0.872678, 0.872678], [-0.872678, 0.872678, 0.872678], [0, -1.41202, 0.539345], [0.872678, -0.872678, 0.872678], [-0.872678, -0.872678, 0.872678], [0.539345, 0, -1.41202], [-0.539345, 0, -1.41202], [0, 1.41202, -0.539345], [0.872678, 0.872678, -0.872678], [-0.872678, 0.872678, -0.872678], [0, -1.41202, -0.539345], [0.872678, -0.872678, -0.872678], [-0.872678, -0.872678, -0.872678], [1.41202, 0.539345, 0], [-1.41202, 0.539345, 0], [1.41202, -0.539345, 0], [-1.41202, -0.539345, 0]]

e = [[1, 2], [1, 4], [1, 7], [2, 5], [2, 8], [3, 4], [3, 5],
  [3, 11], [4, 17], [5, 18], [6, 7], [6, 8], [6, 14], [7, 19],
  [8, 20], [9, 10], [9, 12], [9, 15], [10, 13], [10, 16],
  [11, 12], [11, 13], [12, 17], [13, 18], [14, 15], [14, 16],
  [15, 19], [16, 20], [17, 19], [18, 20]]
Dodecahedron = Polyhedron(v, e, f, atol=1e-4)

####(flexible) candy
f = [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12], [2, 3, 9, 8],
  [3, 4, 10, 9], [5, 6, 12, 11], [1, 6, 12, 7], [2, 8, 13],
  [1, 7, 14], [5, 11, 15], [4, 10, 16], [7, 14, 17],
  [13, 14, 17, 18], [7, 8, 18, 17], [8, 13, 18], [1, 14, 19],
  [13, 14, 19, 20], [1, 2, 20, 19], [2, 13, 20], [10, 16, 21],
  [15, 16, 21, 22], [10, 11, 22, 21], [11, 15, 22], [5, 15, 23],
  [4, 5, 23, 24], [15, 16, 24, 23], [4, 16, 24]]
v = 1.0 * [[0.5, 0.866025, 0.0], [-0.5, 0.866025, 0.0],
  [-1.0, 0.0, 0.0], [-0.5, -0.866025, 0.0], [1 / 2, -0.866025, 0],
  [1.0, 0, 0], [1 / 2, 0.866025, 1], [-1 / 2, 0.866025, 1],
  [-1.0, 0.0, 1], [-1 / 2, -0.866025, 1], [1 / 2, -0.866025, 1],
  [1.0, 0, 1], [-0.5, 1.73205, 0.5], [0.5, 1.73205, 0.5],
  [0.5, -1.73205, 0.5], [-0.5, -1.73205, 0.5],
  [0.5, 1.73205, 1.5], [-0.5, 1.73205, 1.5],
  [0.5, 1.73205, -0.5], [-0.5, 1.73205, -0.5],
  [-0.5, -1.73205, 1.5], [0.5, -1.73205, 1.5],
  [0.5, -1.73205, -0.5], [-0.5, -1.73205, -0.5]]
e = [[1, 2], [1, 6], [1, 7], [1, 14], [1, 19], [2, 3], [2, 8],
  [2, 13], [2, 20], [3, 4], [3, 9], [4, 5], [4, 10], [4, 16],
  [4, 24], [5, 6], [5, 11], [5, 15], [5, 23], [6, 12], [7, 8],
  [7, 12], [7, 14], [7, 17], [8, 9], [8, 13], [8, 18], [9, 10],
  [10, 11], [10, 16], [10, 21], [11, 12], [11, 15], [11, 22],
  [13, 14], [13, 18], [13, 20], [14, 17], [14, 19], [15, 16],
  [15, 22], [15, 23], [16, 21], [16, 24], [17, 18], [19, 20],
  [21, 22], [23, 24]]

Candy = Polyhedron(v, e, f)








# ##################################################################################
# ################# archimedian solids
# ##################################################################################

# ####TruncatedTetrahedron

# TruncatedTetrahedron=Polyhedron(1.0*[ [ 3, -1, -1 ], [ -3, 1, -1 ], [ -3, -1, 1 ], [ 3, 1, 1 ], [ 1, -3, -1 ], [ -1, 3, -1 ], [ -1, -3, 1 ], [ 1, 3, 1 ], [ 1, -1, -3 ], [ -1, 1, -3 ], [ -1, -1, 3 ], [ 1, 1, 3 ] ],[ [ 1, 4 ], [ 1, 5 ], [ 1, 9 ], [ 2, 3 ], [ 2, 6 ], [ 2, 10 ], [ 3, 7 ], [ 3, 11 ], [ 4, 8 ], [ 4, 12 ], [ 5, 7 ], [ 5, 9 ], [ 6, 8 ], [ 6, 10 ], [ 7, 11 ], [ 8, 12 ], [ 9, 10 ], [ 11, 12 ] ],[ [ 1, 5, 9 ], [ 2, 6, 10 ], [ 3, 7, 11 ], [ 4, 8, 12 ], [ 1, 4, 6, 8, 9, 10 ], [ 2, 3, 5, 7, 9, 10 ], [ 2, 3, 6, 8, 11, 12 ], [ 1, 4, 5, 7, 11, 12 ] ])

# ####Cuboctahedron

# Cuboctahedron=Polyhedron([ [ 1, 1, 0 ], [ -1, 1, 0 ], [ 1, -1, 0 ], [ -1, -1, 0 ], [ 1, 0,1 ], [ -1, 0, 1 ], [ 1, 0, -1 ], [ -1, 0, -1 ], [ 0, 1, 1 ], [ 0, -1, 1 ], [ 0, 1, -1 ], [ 0, -1, -1 ] ],[ [ 1, 5 ], [ 1, 7 ], [ 1, 9 ], [ 1, 11 ], [ 2, 6 ], [ 2, 8 ], [ 2, 9 ], [ 2, 11 ], [ 3, 5 ], [ 3, 7 ], [ 3, 10 ], [ 3, 12 ], [ 4, 6 ], [ 4, 8 ], [ 4, 10 ], [ 4, 12 ], [ 5, 9 ], [ 5, 10 ], [ 6, 9 ], [ 6, 10 ], [ 7, 11 ], [ 7, 12 ], [ 8, 11 ], [ 8, 12 ] ],[ [ 1, 5, 9 ], [ 1, 7, 11 ], [ 2, 6, 9 ], [ 2, 8, 11 ], [ 3, 5, 10 ], [ 3, 7, 12 ], [ 4, 6, 10 ], [ 4, 8, 12 ], [ 1, 3, 5, 7 ], [ 5, 6, 9, 10 ], [ 1, 2, 9, 11 ], [ 2, 4, 6, 8 ], [ 3, 4,10, 12 ], [ 7, 8, 11, 12 ] ])

# ####TruncatedCube

# TruncatedCube=Polyhedron([ [ -1., -1., -0.41421400000000003 ], [ -1., -1., 0.41421400000000003 ], [ -1., -0.41421400000000003, -1. ], [ -1., -0.41421400000000003, 1. ], [ -1., 0.41421400000000003, -1. ], [ -1., 0.41421400000000003, 1. ], [ -1., 1., -0.41421400000000003 ], [ -1., 1., 0.41421400000000003 ], [ -0.41421400000000003, -1., -1. ], [ -0.41421400000000003, -1., 1. ], [ -0.41421400000000003, 1., -1. ], [ -0.41421400000000003, 1., 1. ], [ 0.41421400000000003, -1., -1. ], [ 0.41421400000000003, -1., 1. ], [ 0.41421400000000003, 1., -1. ], [ 0.41421400000000003, 1., 1. ], [ 1., -1., -0.41421400000000003 ], [ 1., -1., 0.41421400000000003 ], [ 1., -0.41421400000000003, -1. ], [ 1., -0.41421400000000003,1. ], [ 1., 0.41421400000000003, -1. ], [ 1., 0.41421400000000003, 1. ], [ 1., 1., -0.41421400000000003 ], [ 1., 1., 0.41421400000000003 ], [ 2.7755599999999997e-17, 0., 1. ], [ 1., 2.7755599999999997e-17, 0. ], [ 2.7755599999999997e-17, 0., -1. ], [ -1., 2.7755599999999997e-17, 0. ], [ 2.7755599999999997e-17, 1., 0. ], [ 2.7755599999999997e-17, -1., 0. ] ],[ [ 1, 2 ], [ 1, 3 ], [ 1, 9 ], [ 2, 4 ], [ 2, 10 ], [ 3, 5 ], [ 3, 9 ], [ 4, 6 ], [ 4, 10 ], [ 5, 7 ], [ 5, 11 ], [ 6, 8 ], [ 6, 12 ], [ 7, 8 ], [ 7, 11 ], [ 8, 12 ], [ 9, 13 ], [ 10, 14 ], [ 11, 15 ], [ 12, 16 ], [ 13, 17 ], [ 13, 19 ], [ 14, 18 ], [ 14, 20 ], [ 15, 21 ], [ 15, 23 ], [ 16, 22 ], [ 16, 24 ], [ 17, 18 ], [ 17, 19 ], [ 18, 20 ], [ 19, 21 ], [ 20, 22 ], [ 21, 23 ], [ 22, 24 ], [ 23, 24 ] ],[ [ 1, 3, 9 ], [ 2, 4, 10 ], [ 5, 7, 11 ], [ 6, 8, 12 ], [ 13, 17, 19 ], [ 14, 18, 20 ], [ 15, 21, 23 ], [ 16, 22, 24 ], [ 4, 6, 10, 12, 14, 16, 20, 22 ], [ 17, 18, 19, 20, 21, 22, 23, 24 ], [ 3, 5, 9, 11, 13, 15, 19, 21 ], [ 1, 2, 3, 4, 5, 6, 7, 8 ], [ 7, 8, 11, 12, 15, 16, 23, 24 ], [ 1, 2, 9, 10, 13, 14, 17, 18 ] ])

# ####TruncatedOcathedron
# TruncatedOctahedron=Polyhedron([ [ 2, 1, 0 ], [ 2, 0, 1 ], [ 2, -1, 0 ], [ 2, 0, -1 ], [ -2, 1, 0 ], [ -2, 0, 1 ], [ -2, -1, 0 ], [ -2, 0, -1 ], [ 1, 2, 0 ], [ 0, 2, 1 ], [ -1, 2, 0 ], [ 0, 2, -1 ], [ 1, -2, 0 ], [ 0, -2, 1 ], [ -1, -2, 0 ], [ 0, -2, -1 ], [ 1, 0, 2 ], [ 0, 1, 2 ], [ -1, 0, 2 ], [ 0, -1, 2 ], [ 1, 0, -2 ], [ 0, 1, -2 ], [ -1, 0, -2 ], [ 0, -1, -2 ] ],[ [ 1, 2 ], [ 1, 4 ], [ 1, 9 ], [ 2, 3 ], [ 2, 17 ], [ 3, 4 ], [ 3, 13 ], [ 4, 21 ], [ 5, 6 ], [ 5, 8 ], [ 5, 11 ], [ 6, 7 ], [ 6, 19 ], [ 7, 8 ], [ 7, 15 ], [ 8, 23 ], [ 9, 10 ], [ 9, 12 ], [ 10, 11 ], [ 10, 18 ], [ 11, 12 ], [ 12, 22 ], [ 13, 14 ], [ 13, 16 ], [ 14, 15 ], [ 14, 20 ], [ 15, 16 ], [ 16, 24 ], [ 17, 18 ], [ 17, 20 ], [ 18, 19 ], [ 19, 20 ], [ 21, 22 ], [ 21, 24 ], [ 22, 23 ], [ 23, 24 ] ],[ [ 1, 2, 3, 4 ], [ 5, 6, 7, 8 ], [ 9, 10, 11, 12 ], [ 13, 14, 15, 16 ], [ 17, 18, 19, 20 ], [ 21, 22, 23, 24 ], [ 1, 4, 9, 12, 21, 22 ], [ 3, 4, 13, 16, 21, 24 ], [ 2, 3, 13, 14, 17, 20 ], [ 1, 2, 9, 10, 17, 18 ], [ 5, 6, 10, 11, 18, 19 ], [ 5, 8, 11, 12, 22, 23 ], [ 7, 8, 15, 16, 23, 24 ], [ 6, 7, 14, 15, 19, 20 ] ])

# #Rhombcuboctahedron
# Rhombcuboctahedron=Polyhedron([ [ 1, 1, 2.4142135623730949 ], [ 1, 2.4142135623730949, 1 ], [ 2.4142135623730949, 1, 1 ], [ -1, 1, 2.4142135623730949 ], [ 1, -1, 2.4142135623730949 ], [ 1, 1, -2.4142135623730949 ], [ -1, -1, 2.4142135623730949 ], [ -1, 1, -2.4142135623730949 ], [ 1, -1, -2.4142135623730949 ], [ -1, -1, -2.4142135623730949 ], [ -1, 2.4142135623730949, 1 ], [ 1, -2.4142135623730949, 1 ], [ 1, 2.4142135623730949, -1 ], [ -1, -2.4142135623730949, 1 ], [ -1, 2.4142135623730949, -1 ], [ 1, -2.4142135623730949, -1 ], [ -1, -2.4142135623730949, -1 ], [ -2.4142135623730949, 1, 1 ], [ 2.4142135623730949, -1, 1 ], [ 2.4142135623730949, 1, -1 ], [ -2.4142135623730949, -1, 1 ], [ -2.4142135623730949, 1, -1 ], [ 2.4142135623730949, -1, -1 ], [ -2.4142135623730949, -1, -1 ] ],[ [ 1, 2 ], [ 1, 3 ], [ 1, 4 ], [ 1, 5 ], [ 2, 3 ], [ 2, 11 ], [ 2, 13 ], [ 3, 19 ], [ 3, 20 ], [ 4, 7 ], [ 4, 11 ], [ 4, 18 ], [ 5, 7 ], [ 5, 12 ], [ 5, 19 ], [ 6, 8 ], [ 6, 9 ], [ 6, 13 ], [ 6, 20 ], [ 7, 14 ], [ 7, 21 ], [ 8, 10 ], [ 8, 15 ], [ 8, 22 ], [ 9, 10 ], [ 9, 16 ], [ 9, 23 ], [ 10, 17 ], [ 10, 24 ], [ 11, 15 ], [ 11, 18 ], [ 12, 14 ], [ 12, 16 ], [ 12, 19 ], [ 13, 15 ], [ 13, 20 ], [ 14, 17 ], [ 14, 21 ], [ 15, 22 ], [ 16, 17 ], [ 16, 23 ], [ 17, 24 ], [ 18, 21 ], [ 18, 22 ], [ 19, 23 ], [ 20, 23 ], [ 21, 24 ], [ 22, 24 ] ],[ [ 1, 2, 3 ], [ 4, 11, 18 ], [ 5, 12, 19 ], [ 6, 13, 20 ], [ 7, 14, 21 ], [ 8, 15, 22 ], [ 9, 16, 23 ], [ 10, 17, 24 ], [ 1, 2, 4, 11 ], [ 1, 3, 5, 19 ], [ 1, 4, 5, 7 ], [ 2, 3, 13, 20 ], [ 2, 11, 13, 15 ], [ 3, 19, 20, 23 ], [ 4, 7, 18, 21 ], [ 5, 7, 12, 14 ], [ 6, 8, 9, 10 ], [ 6, 8, 13, 15 ], [ 6, 9, 20, 23 ], [ 8, 10, 22, 24 ], [ 9, 10, 16, 17 ], [ 11, 15, 18, 22 ], [ 12, 14, 16, 17 ], [ 12, 16, 19, 23 ], [ 14, 17, 21, 24 ], [ 18, 21, 22, 24 ] ])

# #### TruncatedCuboctahedron
# TruncatedCuboctahedron=Polyhedron([ [ 1, 2.4142135623730949, 3.8284271247461903 ], [ -1, 2.4142135623730949, 3.8284271247461903 ], [ 1, -2.4142135623730949, 3.8284271247461903 ], [ 1, 2.4142135623730949, -3.8284271247461903 ], [ -1, -2.4142135623730949, 3.8284271247461903 ], [ -1, 2.4142135623730949, -3.8284271247461903 ], [ 1, -2.4142135623730949, -3.8284271247461903 ], [ -1, -2.4142135623730949, -3.8284271247461903 ], [ 3.8284271247461903, 1, 2.4142135623730949 ], [ -3.8284271247461903, 1, 2.4142135623730949 ], [ 3.8284271247461903, -1, 2.4142135623730949 ], [ 3.8284271247461903, 1, -2.4142135623730949 ], [ -3.8284271247461903, -1, 2.4142135623730949 ], [ -3.8284271247461903, 1, -2.4142135623730949 ], [ 3.8284271247461903, -1, -2.4142135623730949 ], [ -3.8284271247461903, -1, -2.4142135623730949 ], [ 2.4142135623730949, 3.8284271247461903, 1 ], [ -2.4142135623730949, 3.8284271247461903, 1 ], [ 2.4142135623730949, -3.8284271247461903, 1 ], [ 2.4142135623730949, 3.8284271247461903, -1 ], [ -2.4142135623730949, -3.8284271247461903, 1 ], [ -2.4142135623730949, 3.8284271247461903, -1 ], [ 2.4142135623730949, -3.8284271247461903, -1 ], [ -2.4142135623730949, -3.8284271247461903, -1 ], [ 2.4142135623730949, 1, 3.8284271247461903 ], [ -2.4142135623730949, 1, 3.8284271247461903 ], [ 2.4142135623730949, -1, 3.8284271247461903 ], [ 2.4142135623730949, 1, -3.8284271247461903 ], [ -2.4142135623730949, -1, 3.8284271247461903 ], [ -2.4142135623730949, 1, -3.8284271247461903 ], [ 2.4142135623730949, -1, -3.8284271247461903 ], [ -2.4142135623730949, -1, -3.8284271247461903 ], [ 1, 3.8284271247461903, 2.4142135623730949 ], [ -1, 3.8284271247461903, 2.4142135623730949 ], [ 1, -3.8284271247461903, 2.4142135623730949 ], [ 1, 3.8284271247461903, -2.4142135623730949 ], [ -1, -3.8284271247461903, 2.4142135623730949 ], [ -1, 3.8284271247461903, -2.4142135623730949 ], [ 1, -3.8284271247461903, -2.4142135623730949 ], [ -1, -3.8284271247461903, -2.4142135623730949 ], [ 3.8284271247461903, 2.4142135623730949, 1 ], [ -3.8284271247461903, 2.4142135623730949, 1 ], [ 3.8284271247461903, -2.4142135623730949, 1 ], [ 3.8284271247461903, 2.4142135623730949, -1 ], [ -3.8284271247461903, -2.4142135623730949, 1 ], [ -3.8284271247461903, 2.4142135623730949, -1 ], [ 3.8284271247461903, -2.4142135623730949, -1 ], [ -3.8284271247461903, -2.4142135623730949, -1 ] ],[ [ 1, 2 ], [ 1, 25 ], [ 1, 33 ], [ 2, 26 ], [ 2, 34 ], [ 3, 5 ], [ 3, 27 ], [ 3, 35 ], [ 4, 6 ], [ 4, 28 ], [ 4, 36 ], [ 5, 29 ], [ 5, 37 ], [ 6, 30 ], [ 6, 38 ], [ 7, 8 ], [ 7, 31 ], [ 7, 39 ], [ 8, 32 ], [ 8, 40 ], [ 9, 11 ], [ 9, 25 ], [ 9, 41 ], [ 10, 13 ], [ 10, 26 ], [ 10, 42 ], [ 11, 27 ], [ 11, 43 ], [ 12, 15 ], [ 12, 28 ], [ 12, 44 ], [ 13, 29 ], [ 13, 45 ], [ 14, 16 ], [ 14, 30 ], [ 14, 46 ], [ 15, 31 ], [ 15, 47 ], [ 16, 32 ], [ 16, 48 ], [ 17, 20 ], [ 17, 33 ], [ 17, 41 ], [ 18, 22 ], [ 18, 34 ], [ 18, 42 ], [ 19, 23 ], [ 19, 35 ], [ 19, 43 ], [ 20, 36 ], [ 20, 44 ], [ 21, 24 ], [ 21, 37 ], [ 21, 45 ], [ 22, 38 ], [ 22, 46 ], [ 23, 39 ], [ 23, 47 ], [ 24, 40 ], [ 24, 48 ], [ 25, 27 ], [ 26, 29 ], [ 28, 31 ], [ 30, 32 ], [ 33, 34 ], [ 35, 37 ], [ 36, 38 ], [ 39, 40 ], [ 41, 44 ], [ 42, 46 ], [ 43, 47 ], [ 45, 48 ] ],[ [ 1, 2, 33, 34 ], [ 3, 5, 35, 37 ], [ 4, 6, 36, 38 ], [ 7, 8, 39, 40 ], [ 9, 11, 25, 27 ], [ 10, 13, 26, 29 ], [ 12, 15, 28, 31 ], [ 14, 16, 30, 32 ], [ 17, 20, 41, 44 ], [ 18, 22, 42, 46 ], [ 19, 23, 43, 47 ], [ 21, 24, 45, 48 ], [ 1, 9, 17, 25, 33, 41 ], [ 2, 10, 18, 26, 34, 42 ], [ 3, 11, 19, 27, 35, 43 ], [ 4, 12, 20, 28, 36, 44 ], [ 5, 13, 21, 29, 37, 45 ], [ 6, 14, 22, 30, 38, 46 ], [ 7, 15, 23, 31, 39, 47 ], [ 8, 16, 24, 32, 40, 48 ], [ 1, 2, 3, 5, 25, 26, 27, 29 ], [ 19, 21, 23, 24, 35, 37, 39, 40 ], [ 4, 6, 7, 8, 28, 30, 31, 32 ], [ 17, 18, 20, 22, 33, 34, 36, 38 ], [ 9, 11, 12, 15, 41, 43, 44, 47 ], [ 10, 13, 14, 16, 42, 45, 46, 48 ] ])

# ####Snubcube
# Snubcube=Polyhedron([ [ 0.54368901269207637, 1., 1.839286755214161 ], [ 1.839286755214161, 0.54368901269207637, 1. ], [ 1., 1.839286755214161, 0.54368901269207637 ], [ -1., 0.54368901269207637, 1.839286755214161 ], [ 1.839286755214161, -1., 0.54368901269207637 ], [ 0.54368901269207637, 1.839286755214161, -1. ], [ 1., -0.54368901269207637, 1.839286755214161 ], [ 1.839286755214161, 1., -0.54368901269207637 ], [ -0.54368901269207637, 1.839286755214161, 1. ], [ 1., 0.54368901269207637, -1.839286755214161 ], [ -1.839286755214161, 1., 0.54368901269207637 ], [ 0.54368901269207637, -1.839286755214161, 1. ], [ -0.54368901269207637, -1., 1.839286755214161 ], [ 1.839286755214161, -0.54368901269207637, -1. ], [ -1., 1.839286755214161, -0.54368901269207637 ], [ 0.54368901269207637, -1., -1.839286755214161 ], [ -1.839286755214161, 0.54368901269207637, -1. ], [ -1., -1.839286755214161, 0.54368901269207637 ], [ -0.54368901269207637, 1., -1.839286755214161 ], [ -1.839286755214161, -0.54368901269207637, 1. ], [ 1., -1.839286755214161, -0.54368901269207637 ], [ -1., -0.54368901269207637, -1.839286755214161 ], [ -1.839286755214161, -1., -0.54368901269207637 ], [ -0.54368901269207637, -1.839286755214161, -1. ] ],[ [ 1, 2 ], [ 1, 3 ], [ 1, 4 ], [ 1, 7 ], [ 1, 9 ], [ 2, 3 ], [ 2, 5 ], [ 2, 7 ], [ 2, 8 ], [ 3, 6 ], [ 3, 8 ], [ 3, 9 ], [ 4, 9 ], [ 4, 11 ], [ 4, 13 ], [ 4, 20 ], [ 5, 7 ], [ 5, 12 ], [ 5, 14 ], [ 5, 21 ], [ 6, 8 ], [ 6, 10 ], [ 6, 15 ], [ 6, 19 ], [ 7, 12 ], [ 7, 13 ], [ 8, 10 ], [ 8, 14 ], [ 9, 11 ], [ 9, 15 ], [ 10, 14 ], [ 10, 16 ], [ 10, 19 ], [ 11, 15 ], [ 11, 17 ], [ 11, 20 ], [ 12, 13 ], [ 12, 18 ], [ 12, 21 ], [ 13, 18 ], [ 13, 20 ], [ 14, 16 ], [ 14, 21 ], [ 15, 17 ], [ 15, 19 ], [ 16, 21 ], [ 16, 22 ], [ 16, 24 ], [ 17, 19 ], [ 17, 22 ], [ 17, 23 ], [ 18, 20 ], [ 18, 23 ], [ 18, 24 ], [ 19, 22 ], [ 20, 23 ], [ 21, 24 ], [ 22, 23 ], [ 22, 24 ], [ 23, 24 ] ],[ [ 1, 2, 3 ], [ 1, 2, 7 ], [ 1, 3, 9 ], [ 1, 4, 9 ], [ 2, 3, 8 ], [ 2, 5, 7 ], [ 3, 6, 8 ], [ 4, 9, 11 ], [ 4, 11, 20 ], [ 4, 13, 20 ], [ 5, 7, 12 ], [ 5, 12, 21 ], [ 5, 14, 21 ], [ 6, 8, 10 ], [ 6, 10, 19 ], [ 6, 15, 19 ], [ 7, 12, 13 ], [ 8, 10, 14 ], [ 9, 11, 15 ], [ 10, 14, 16 ], [ 11, 15, 17 ], [ 12, 13, 18 ], [ 13, 18, 20 ], [ 14, 16, 21 ], [ 15, 17, 19 ], [ 16, 21, 24 ], [ 16, 22, 24 ], [ 17, 19, 22 ], [ 17, 22, 23 ], [ 18, 20, 23 ], [ 18, 23, 24 ], [ 22, 23, 24 ], [ 11, 17, 20, 23 ], [ 10, 16, 19, 22 ], [ 12, 18, 21, 24 ], [ 1, 4, 7, 13 ], [ 3, 6, 9, 15 ], [ 2, 5, 8, 14 ] ])

# ####Icosidodecahedron
#  Icosidodecahedron=Polyhedron([ [ 0., 0., 1.6180339887498949 ], [ 0., 1.6180339887498949, 0. ], [ 1.6180339887498949, 0., 0. ], [ 0., 0., -1.6180339887498949 ], [ 0., -1.6180339887498949, 0. ], [ -1.6180339887498949, 0., 0. ], [ 0.5, 0.80901699437494745, 1.3090169943749475 ], [ 0.80901699437494745, 1.3090169943749475, 0.5 ], [ 1.3090169943749475, 0.5, 0.80901699437494745 ], [ -0.5, 0.80901699437494745, 1.3090169943749475 ], [ 0.80901699437494745, 1.3090169943749475, -0.5 ], [ 1.3090169943749475, -0.5, 0.80901699437494745 ], [ 0.5, -0.80901699437494745, 1.3090169943749475 ], [ -0.80901699437494745, 1.3090169943749475, 0.5 ], [ 1.3090169943749475, 0.5, -0.80901699437494745 ], [ 0.5, 0.80901699437494745, -1.3090169943749475 ], [ 0.80901699437494745, -1.3090169943749475, 0.5 ], [ -1.3090169943749475, 0.5, 0.80901699437494745 ], [ -0.5, -0.80901699437494745, 1.3090169943749475 ], [ -0.80901699437494745, 1.3090169943749475, -0.5 ], [ 1.3090169943749475, -0.5, -0.80901699437494745 ], [ 0.5, -0.80901699437494745, -1.3090169943749475 ], [ -0.80901699437494745, -1.3090169943749475, 0.5 ], [ -1.3090169943749475, 0.5, -0.80901699437494745 ], [ -0.5, 0.80901699437494745, -1.3090169943749475 ], [ 0.80901699437494745, -1.3090169943749475, -0.5 ], [ -1.3090169943749475, -0.5, 0.80901699437494745 ], [ -0.5, -0.80901699437494745, -1.3090169943749475 ], [ -0.80901699437494745, -1.3090169943749475, -0.5 ], [ -1.3090169943749475, -0.5, -0.80901699437494745 ] ],[ [ 1, 7 ], [ 1, 10 ], [ 1, 13 ], [ 1, 19 ], [ 2, 8 ], [ 2, 11 ], [ 2, 14 ], [ 2, 20 ], [ 3, 9 ], [ 3, 12 ], [ 3, 15 ], [ 3, 21 ], [ 4, 16 ], [ 4, 22 ], [ 4, 25 ], [ 4, 28 ], [ 5, 17 ], [ 5, 23 ], [ 5, 26 ], [ 5, 29 ], [ 6, 18 ], [ 6, 24 ], [ 6, 27 ], [ 6, 30 ], [ 7, 8 ], [ 7, 9 ], [ 7, 10 ], [ 8, 9 ], [ 8, 11 ], [ 9, 12 ], [ 10, 14 ], [ 10, 18 ], [ 11, 15 ], [ 11, 16 ], [ 12, 13 ], [ 12, 17 ], [ 13, 17 ], [ 13, 19 ], [ 14, 18 ], [ 14, 20 ], [ 15, 16 ], [ 15, 21 ], [ 16, 25 ], [ 17, 26 ], [ 18, 27 ], [ 19, 23 ], [ 19, 27 ], [ 20, 24 ], [ 20, 25 ], [ 21, 22 ], [ 21, 26 ], [ 22, 26 ], [ 22, 28 ], [ 23, 27 ], [ 23, 29 ], [ 24, 25 ], [ 24, 30 ], [ 28, 29 ], [ 28, 30 ], [ 29, 30 ] ],[ [ 1, 7, 10 ], [ 1, 13, 19 ], [ 2, 8, 11 ], [ 2, 14, 20 ], [ 3, 9, 12 ], [ 3, 15, 21 ], [ 4, 16, 25 ], [ 4, 22, 28 ], [ 5, 17, 26 ], [ 5, 23, 29 ], [ 6, 18, 27 ], [ 6, 24, 30 ], [ 7, 8, 9 ], [ 10, 14, 18 ], [ 11, 15, 16 ], [ 12, 13, 17 ], [ 19, 23, 27 ], [ 20, 24, 25 ], [ 21, 22, 26 ], [ 28, 29, 30 ], [ 1, 7, 9, 12, 13 ], [ 1, 10, 18, 19, 27 ], [ 2, 7, 8, 10, 14 ], [ 2, 11, 16, 20, 25 ], [ 3, 8, 9, 11, 15 ], [ 3, 12, 17, 21, 26 ], [ 4, 15, 16, 21, 22 ], [ 4, 24, 25, 28, 30 ], [ 5, 13, 17, 19, 23 ], [ 5, 22, 26, 28, 29 ], [ 6, 14, 18, 20, 24 ], [ 6, 23, 27, 29, 30 ] ])

# ####TruncatedDodecahedron
#  TruncatedDodecahedron=Polyhedron([ [ 0, 0.61803398874989479, 3.6180339887498949 ], [ 3.6180339887498949, 0, 0.61803398874989479 ], [ 0.61803398874989479, 3.6180339887498949, 0 ], [ 0, -0.61803398874989479, 3.6180339887498949 ], [ 3.6180339887498949, 0, -0.61803398874989479 ], [ -0.61803398874989479, 3.6180339887498949, 0 ], [ 0, 0.61803398874989479, -3.6180339887498949 ], [ -3.6180339887498949, 0, 0.61803398874989479 ], [ 0.61803398874989479, -3.6180339887498949, 0 ], [ 0, -0.61803398874989479, -3.6180339887498949 ], [ -3.6180339887498949, 0, -0.61803398874989479 ], [ -0.61803398874989479, -3.6180339887498949, 0 ], [ 0.61803398874989479, 1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, 0.61803398874989479, 1.6180339887498949 ], [ 1.6180339887498949, 3.2360679774997898, 0.61803398874989479 ], [ -0.61803398874989479, 1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, -0.61803398874989479, 1.6180339887498949 ], [ 1.6180339887498949, 3.2360679774997898, -0.61803398874989479 ], [ 0.61803398874989479, -1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, 0.61803398874989479, -1.6180339887498949 ], [ -1.6180339887498949, 3.2360679774997898, 0.61803398874989479 ], [ 0.61803398874989479, 1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, 0.61803398874989479, 1.6180339887498949 ], [ 1.6180339887498949, -3.2360679774997898, 0.61803398874989479 ], [ -0.61803398874989479, -1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, -0.61803398874989479, -1.6180339887498949 ], [ -1.6180339887498949, 3.2360679774997898, -0.61803398874989479 ], [ -0.61803398874989479, 1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, -0.61803398874989479, 1.6180339887498949 ], [ 1.6180339887498949, -3.2360679774997898, -0.61803398874989479 ], [ 0.61803398874989479, -1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, 0.61803398874989479, -1.6180339887498949 ], [ -1.6180339887498949, -3.2360679774997898, 0.61803398874989479 ], [ -0.61803398874989479, -1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, -0.61803398874989479, -1.6180339887498949 ], [ -1.6180339887498949, -3.2360679774997898, -0.61803398874989479 ], [ 1.6180339887498949, 2, 2.6180339887498949 ], [ 2.6180339887498949, 1.6180339887498949, 2 ], [ 2, 2.6180339887498949, 1.6180339887498949 ], [ -1.6180339887498949, 2, 2.6180339887498949 ], [ 2.6180339887498949, -1.6180339887498949, 2 ], [ 2, 2.6180339887498949, -1.6180339887498949 ], [ 1.6180339887498949, -2, 2.6180339887498949 ], [ 2.6180339887498949, 1.6180339887498949, -2 ], [ -2, 2.6180339887498949, 1.6180339887498949 ], [ 1.6180339887498949, 2, -2.6180339887498949 ], [ -2.6180339887498949, 1.6180339887498949, 2 ], [ 2, -2.6180339887498949, 1.6180339887498949 ], [ -1.6180339887498949, -2, 2.6180339887498949 ], [ 2.6180339887498949, -1.6180339887498949, -2 ], [ -2, 2.6180339887498949, -1.6180339887498949 ], [ -1.6180339887498949, 2, -2.6180339887498949 ], [ -2.6180339887498949, -1.6180339887498949, 2 ], [ 2, -2.6180339887498949, -1.6180339887498949 ], [ 1.6180339887498949, -2, -2.6180339887498949 ], [ -2.6180339887498949, 1.6180339887498949, -2 ], [ -2, -2.6180339887498949, 1.6180339887498949 ], [ -1.6180339887498949, -2, -2.6180339887498949 ], [ -2.6180339887498949, -1.6180339887498949, -2 ], [ -2, -2.6180339887498949, -1.6180339887498949 ] ],[ [ 1, 4 ], [ 1, 13 ], [ 1, 16 ], [ 2, 5 ], [ 2, 14 ], [ 2, 17 ], [ 3, 6 ], [ 3, 15 ], [ 3, 18 ], [ 4, 19 ], [ 4, 25 ], [ 5, 20 ], [ 5, 26 ], [ 6, 21 ], [ 6, 27 ], [ 7, 10 ], [ 7, 22 ], [ 7, 28 ], [ 8, 11 ], [ 8, 23 ], [ 8, 29 ], [ 9, 12 ], [ 9, 24 ], [ 9, 30 ], [ 10, 31 ], [ 10, 34 ], [ 11, 32 ], [ 11, 35 ], [ 12, 33 ], [ 12, 36 ], [ 13, 16 ], [ 13, 37 ], [ 14, 17 ], [ 14, 38 ], [ 15, 18 ], [ 15, 39 ], [ 16, 40 ], [ 17, 41 ], [ 18, 42 ], [ 19, 25 ], [ 19, 43 ], [ 20, 26 ], [ 20, 44 ], [ 21, 27 ], [ 21, 45 ], [ 22, 28 ], [ 22, 46 ], [ 23, 29 ], [ 23, 47 ], [ 24, 30 ], [ 24, 48 ], [ 25, 49 ], [ 26, 50 ], [ 27, 51 ], [ 28, 52 ], [ 29, 53 ], [ 30, 54 ], [ 31, 34 ], [ 31, 55 ], [ 32, 35 ], [ 32, 56 ], [ 33, 36 ], [ 33, 57 ], [ 34, 58 ], [ 35, 59 ], [ 36, 60 ], [ 37, 38 ], [ 37, 39 ], [ 38, 39 ], [ 40, 45 ], [ 40, 47 ], [ 41, 43 ], [ 41, 48 ], [ 42, 44 ], [ 42, 46 ], [ 43, 48 ], [ 44, 46 ], [ 45, 47 ], [ 49, 53 ], [ 49, 57 ], [ 50, 54 ], [ 50, 55 ], [ 51, 52 ], [ 51, 56 ], [ 52, 56 ], [ 53, 57 ], [ 54, 55 ], [ 58, 59 ], [ 58, 60 ], [ 59, 60 ] ],[ [ 1, 13, 16 ], [ 2, 14, 17 ], [ 3, 15, 18 ], [ 4, 19, 25 ], [ 5, 20, 26 ], [ 6, 21, 27 ], [ 7, 22, 28 ], [ 8, 23, 29 ], [ 9, 24, 30 ], [ 10, 31, 34 ], [ 11, 32, 35 ], [ 12, 33, 36 ], [ 37, 38, 39 ], [ 40, 45, 47 ], [ 41, 43, 48 ], [ 42, 44, 46 ], [ 49, 53, 57 ], [ 50, 54, 55 ], [ 51, 52, 56 ], [ 58, 59, 60 ], [ 1, 4, 13, 14, 17, 19, 37, 38, 41, 43 ], [ 1, 4, 16, 23, 25, 29, 40, 47, 49, 53 ], [ 2, 5, 14, 15, 18, 20, 38, 39, 42, 44 ], [ 2, 5, 17, 24, 26, 30, 41, 48, 50, 54 ], [ 3, 6, 13, 15, 16, 21, 37, 39, 40, 45 ], [ 3, 6, 18, 22, 27, 28, 42, 46, 51, 52 ], [ 7, 10, 20, 22, 26, 31, 44, 46, 50, 55 ], [ 7, 10, 28, 32, 34, 35, 52, 56, 58, 59 ], [ 8, 11, 21, 23, 27, 32, 45, 47, 51, 56 ], [ 8, 11, 29, 33, 35, 36, 53, 57, 59, 60 ], [ 9, 12, 19, 24, 25, 33, 43, 48, 49, 57 ], [ 9, 12, 30, 31, 34, 36, 54, 55, 58, 60 ] ])

# ####TruncatedIcosahedron
#  TruncatedIcosahedron=Polyhedron([ [ 0, 1, 4.8541019662496847 ], [ 4.8541019662496847, 0, 1 ], [ 1, 4.8541019662496847, 0 ], [ 0, -1, 4.8541019662496847 ], [ 4.8541019662496847, 0, -1 ], [ -1, 4.8541019662496847, 0 ], [ 0, 1, -4.8541019662496847 ], [ -4.8541019662496847, 0, 1 ], [ 1, -4.8541019662496847, 0 ], [ 0, -1, -4.8541019662496847 ], [ -4.8541019662496847, 0, -1 ], [ -1, -4.8541019662496847, 0 ], [ 1, 3.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, 1, 3.6180339887498949 ], [ 3.6180339887498949, 3.2360679774997898, 1 ], [ -1, 3.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, -1, 3.6180339887498949 ], [ 3.6180339887498949, 3.2360679774997898, -1 ], [ 1, -3.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, 1, -3.6180339887498949 ], [ -3.6180339887498949, 3.2360679774997898, 1 ], [ 1, 3.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, 1, 3.6180339887498949 ], [ 3.6180339887498949, -3.2360679774997898, 1 ], [ -1, -3.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, -1, -3.6180339887498949 ], [ -3.6180339887498949, 3.2360679774997898, -1 ], [ -1, 3.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, -1, 3.6180339887498949 ], [ 3.6180339887498949, -3.2360679774997898, -1 ], [ 1, -3.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, 1, -3.6180339887498949 ], [ -3.6180339887498949, -3.2360679774997898, 1 ], [ -1, -3.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, -1, -3.6180339887498949 ], [ -3.6180339887498949, -3.2360679774997898, -1 ], [ 1.6180339887498949, 2, 4.2360679774997898 ], [ 4.2360679774997898, 1.6180339887498949, 2 ], [ 2, 4.2360679774997898, 1.6180339887498949 ], [ -1.6180339887498949, 2, 4.2360679774997898 ], [ 4.2360679774997898, -1.6180339887498949, 2 ], [ 2, 4.2360679774997898, -1.6180339887498949 ], [ 1.6180339887498949, -2, 4.2360679774997898 ], [ 4.2360679774997898, 1.6180339887498949, -2 ], [ -2, 4.2360679774997898, 1.6180339887498949 ], [ 1.6180339887498949, 2, -4.2360679774997898 ], [ -4.2360679774997898, 1.6180339887498949, 2 ], [ 2, -4.2360679774997898, 1.6180339887498949 ], [ -1.6180339887498949, -2, 4.2360679774997898 ], [ 4.2360679774997898, -1.6180339887498949, -2 ], [ -2, 4.2360679774997898, -1.6180339887498949 ], [ -1.6180339887498949, 2, -4.2360679774997898 ], [ -4.2360679774997898, -1.6180339887498949, 2 ], [ 2, -4.2360679774997898, -1.6180339887498949 ], [ 1.6180339887498949, -2, -4.2360679774997898 ], [ -4.2360679774997898, 1.6180339887498949, -2 ], [ -2, -4.2360679774997898, 1.6180339887498949 ], [ -1.6180339887498949, -2, -4.2360679774997898 ], [ -4.2360679774997898, -1.6180339887498949, -2 ], [ -2, -4.2360679774997898, -1.6180339887498949 ] ],[ [ 1, 4 ], [ 1, 37 ], [ 1, 40 ], [ 2, 5 ], [ 2, 38 ], [ 2, 41 ], [ 3, 6 ], [ 3, 39 ], [ 3, 42 ], [ 4, 43 ], [ 4, 49 ], [ 5, 44 ], [ 5, 50 ], [ 6, 45 ], [ 6, 51 ], [ 7, 10 ], [ 7, 46 ], [ 7, 52 ], [ 8, 11 ], [ 8, 47 ], [ 8, 53 ], [ 9, 12 ], [ 9, 48 ], [ 9, 54 ], [ 10, 55 ], [ 10, 58 ], [ 11, 56 ], [ 11, 59 ], [ 12, 57 ], [ 12, 60 ], [ 13, 16 ], [ 13, 37 ], [ 13, 39 ], [ 14, 17 ], [ 14, 37 ], [ 14, 38 ], [ 15, 18 ], [ 15, 38 ], [ 15, 39 ], [ 16, 40 ], [ 16, 45 ], [ 17, 41 ], [ 17, 43 ], [ 18, 42 ], [ 18, 44 ], [ 19, 25 ], [ 19, 43 ], [ 19, 48 ], [ 20, 26 ], [ 20, 44 ], [ 20, 46 ], [ 21, 27 ], [ 21, 45 ], [ 21, 47 ], [ 22, 28 ], [ 22, 42 ], [ 22, 46 ], [ 23, 29 ], [ 23, 40 ], [ 23, 47 ], [ 24, 30 ], [ 24, 41 ], [ 24, 48 ], [ 25, 49 ], [ 25, 57 ], [ 26, 50 ], [ 26, 55 ], [ 27, 51 ], [ 27, 56 ], [ 28, 51 ], [ 28, 52 ], [ 29, 49 ], [ 29, 53 ], [ 30, 50 ], [ 30, 54 ], [ 31, 34 ], [ 31, 54 ], [ 31, 55 ], [ 32, 35 ], [ 32, 52 ], [ 32, 56 ], [ 33, 36 ], [ 33, 53 ], [ 33, 57 ], [ 34, 58 ], [ 34, 60 ], [ 35, 58 ], [ 35, 59 ], [ 36, 59 ], [ 36, 60 ] ],[ [ 1, 13, 16, 37, 40 ], [ 2, 14, 17, 38, 41 ], [ 3, 15, 18, 39, 42 ], [ 4, 19, 25, 43, 49 ], [ 5, 20, 26, 44, 50 ], [ 6, 21, 27, 45, 51 ], [ 7, 22, 28, 46, 52 ], [ 8, 23, 29, 47, 53 ], [ 9, 24, 30, 48, 54 ], [ 10, 31, 34, 55, 58 ], [ 11, 32, 35, 56, 59 ], [ 12, 33, 36, 57, 60 ], [ 1, 4, 14, 17, 37, 43 ], [ 1, 4, 23, 29, 40, 49 ], [ 2, 5, 15, 18, 38, 44 ], [ 2, 5, 24, 30, 41, 50 ], [ 3, 6, 13, 16, 39, 45 ], [ 3, 6, 22, 28, 42, 51 ], [ 7, 10, 20, 26, 46, 55 ], [ 7, 10, 32, 35, 52, 58 ], [ 8, 11, 21, 27, 47, 56 ], [ 8, 11, 33, 36, 53, 59 ], [ 9, 12, 19, 25, 48, 57 ], [ 9, 12, 31, 34, 54, 60 ], [ 13, 14, 15, 37, 38, 39 ], [ 16, 21, 23, 40, 45, 47 ], [ 17, 19, 24, 41, 43, 48 ], [ 18, 20, 22, 42, 44, 46 ], [ 25, 29, 33, 49, 53, 57 ], [ 26, 30, 31, 50, 54, 55 ], [ 27, 28, 32, 51, 52, 56 ], [ 34, 35, 36, 58, 59, 60 ] ])

# ####Rhombicosidodecahedron
#  Rhombicosidodecahedron=Polyhedron([ [ 1, 1, 4.2360679774997898 ], [ 4.2360679774997898, 1, 1 ], [ 1, 4.2360679774997898, 1 ], [ -1, 1, 4.2360679774997898 ], [ 4.2360679774997898, -1, 1 ], [ 1, 4.2360679774997898, -1 ], [ 1, -1, 4.2360679774997898 ], [ 4.2360679774997898, 1, -1 ], [ -1, 4.2360679774997898, 1 ], [ 1, 1, -4.2360679774997898 ], [ -4.2360679774997898, 1, 1 ], [ 1, -4.2360679774997898, 1 ], [ -1, -1, 4.2360679774997898 ], [ 4.2360679774997898, -1, -1 ], [ -1, 4.2360679774997898, -1 ], [ -1, 1, -4.2360679774997898 ], [ -4.2360679774997898, -1, 1 ], [ 1, -4.2360679774997898, -1 ], [ 1, -1, -4.2360679774997898 ], [ -4.2360679774997898, 1, -1 ], [ -1, -4.2360679774997898, 1 ], [ -1, -1, -4.2360679774997898 ], [ -4.2360679774997898, -1, -1 ], [ -1, -4.2360679774997898, -1 ], [ 2.6180339887498949, 1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, 2.6180339887498949, 1.6180339887498949 ], [ 1.6180339887498949, 3.2360679774997898, 2.6180339887498949 ], [ -2.6180339887498949, 1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, -2.6180339887498949, 1.6180339887498949 ], [ 1.6180339887498949, 3.2360679774997898, -2.6180339887498949 ], [ 2.6180339887498949, -1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, 2.6180339887498949, -1.6180339887498949 ], [ -1.6180339887498949, 3.2360679774997898, 2.6180339887498949 ], [ 2.6180339887498949, 1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, 2.6180339887498949, 1.6180339887498949 ], [ 1.6180339887498949, -3.2360679774997898, 2.6180339887498949 ], [ -2.6180339887498949, -1.6180339887498949, 3.2360679774997898 ], [ 3.2360679774997898, -2.6180339887498949, -1.6180339887498949 ], [ -1.6180339887498949, 3.2360679774997898, -2.6180339887498949 ], [ -2.6180339887498949, 1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, -2.6180339887498949, 1.6180339887498949 ], [ 1.6180339887498949, -3.2360679774997898, -2.6180339887498949 ], [ 2.6180339887498949, -1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, 2.6180339887498949, -1.6180339887498949 ], [ -1.6180339887498949, -3.2360679774997898, 2.6180339887498949 ], [ -2.6180339887498949, -1.6180339887498949, -3.2360679774997898 ], [ -3.2360679774997898, -2.6180339887498949, -1.6180339887498949 ], [ -1.6180339887498949, -3.2360679774997898, -2.6180339887498949 ], [ 3.6180339887498949, 0, 2.6180339887498949 ], [ 2.6180339887498949, 3.6180339887498949, 0 ], [ 0, 2.6180339887498949, 3.6180339887498949 ], [ -3.618033988498949, 0, 2.6180339887498949 ], [ 2.6180339887498949, -3.6180339887498949, 0 ], [ 0, 2.6180339887498949, -3.6180339887498949 ], [ 3.6180339887498949, 0, -2.6180339887498949 ], [ -2.6180339887498949, 3.6180339887498949, 0 ], [ 0, -2.6180339887498949, 3.6180339887498949 ], [ -3.6180339887498949, 0, -2.6180339887498949 ], [ -2.6180339887498949, -3.6180339887498949, 0 ], [ 0, -2.6180339887498949, -3.6180339887498949 ] ],[ [ 1, 4 ], [ 1, 7 ], [ 1, 25 ], [ 1, 51 ], [ 2, 5 ], [ 2, 8 ], [ 2, 26 ], [ 2, 49 ], [ 3, 6 ], [ 3, 9 ], [ 3, 27 ], [ 3, 50 ], [ 4, 13 ], [ 4, 28 ], [ 4, 51 ], [ 5, 14 ], [ 5, 29 ], [ 5, 49 ], [ 6, 15 ], [ 6, 30 ], [ 6, 50 ], [ 7, 13 ], [ 7, 31 ], [ 7, 57 ], [ 8, 14 ], [ 8, 32 ], [ 8, 55 ], [ 9, 15 ], [ 9, 33 ], [ 9, 56 ], [ 10, 16 ], [ 10, 19 ], [ 10, 34 ], [ 10, 54 ], [ 11, 17 ], [ 11, 20 ], [ 11, 35 ], [ 11, 52 ], [ 12, 18 ], [12, 21 ], [ 12, 36 ], [ 12, 53 ], [ 13, 37 ], [ 13, 57 ], [ 14, 38 ], [ 14, 55 ], [ 15, 39 ], [ 15, 56 ], [ 16, 22 ], [ 16, 40 ], [ 16, 54 ], [ 17, 23 ], [ 17, 41 ], [ 17, 52 ], [ 18, 24 ], [ 18, 42 ], [ 18, 53 ], [ 19, 22 ], [ 19, 43 ], [ 19, 60 ], [ 20, 23 ], [ 20, 44 ], [ 20, 58 ], [ 21, 24 ], [ 21, 45 ], [ 21, 59 ], [ 22, 46 ], [ 22, 60 ], [ 23, 47 ], [ 23, 58 ], [ 24, 48 ], [ 24, 59 ], [ 25, 26 ], [ 25, 27 ], [ 25, 49 ], [ 26, 27 ], [ 26, 50 ], [ 27, 51 ], [ 28, 33 ], [ 28, 35 ], [ 28, 52 ], [ 29, 31 ], [ 29, 36 ], [ 29, 53 ], [ 30, 32 ], [ 30, 34 ], [ 30, 54 ], [ 31, 36 ], [ 31, 49 ], [ 32, 34 ], [ 32, 50 ], [ 33, 35 ], [ 33, 51 ], [ 34, 55 ], [ 35, 56 ], [ 36, 57 ], [ 37, 41 ], [ 37, 45 ], [ 37, 52 ], [ 38, 42 ], [ 38, 43 ], [ 38, 53 ], [ 39, 40 ], [ 39, 44 ], [ 39, 54 ], [ 40, 44 ], [ 40, 58 ], [ 41, 45 ], [ 41, 59 ], [ 42, 43 ], [ 42, 60 ], [ 43, 55 ], [ 44, 56 ], [ 45, 57 ], [ 46, 47 ], [ 46, 48 ], [ 46, 58 ], [ 47, 48 ], [ 47, 59 ], [ 48, 60 ] ],[ [ 1, 4, 51 ], [ 2, 5, 49 ], [ 3, 6, 50 ], [ 7, 13, 57 ], [ 8, 14, 55 ], [ 9, 15, 56 ], [ 10, 16, 54 ], [ 11, 17, 52 ], [ 12, 18, 53 ], [ 19, 22, 60 ], [ 20, 23, 58 ], [ 21, 24, 59 ], [ 25, 26, 27 ], [ 28, 33, 35 ], [ 29, 31, 36 ], [ 30, 32, 34 ], [ 37, 41, 45 ], [ 38, 42, 43 ], [ 39, 40, 44 ], [ 46, 47, 48 ], [ 1, 4, 7, 13 ], [ 1, 25, 27, 51 ], [ 2, 5, 8, 14 ], [ 2, 25, 26, 49 ], [ 3, 6, 9, 15 ], [ 3, 26, 27, 50 ], [ 4, 28, 33, 51 ], [ 5, 29, 31, 49 ], [ 6, 30, 32, 50 ], [ 7, 31, 36, 57 ], [ 8, 32, 34, 55 ], [ 9, 33, 35, 56 ], [ 10, 16, 19, 22 ], [ 10, 30, 34, 54 ], [ 11, 17, 20, 23 ], [ 11, 28, 35, 52 ], [ 12, 18, 21, 24 ], [ 12, 29, 36, 53 ], [ 13, 37, 45, 57 ], [ 14, 38, 43, 55 ], [ 15, 39, 44, 56 ], [ 16, 39, 40, 54 ], [ 17, 37, 41, 52 ], [ 18, 38, 42, 53 ], [ 19, 42, 43, 60 ], [ 20, 40, 44, 58 ], [ 21, 41, 45, 59 ], [ 22, 46, 48, 60 ], [ 23, 46, 47, 58 ], [ 24, 47, 48, 59 ], [ 1, 7, 25, 31, 49 ], [ 2, 8, 26, 32, 50 ], [ 3, 9, 27, 33, 51 ], [ 4, 13, 28, 37, 52 ], [ 5, 14, 29, 38, 53 ], [ 6, 15, 30, 39, 54 ], [ 10, 19, 34, 43, 55 ], [ 11, 20, 35, 44, 56 ], [ 12, 21, 36, 45, 57 ], [ 16, 22, 40, 46, 58 ], [ 17, 23, 41, 47, 59 ], [ 18, 24, 42, 48, 60 ] ])

# ####Truncatedicosidodecahedron
#  Truncatedicosidodecahedron=Polyhedron([ [ 0.61803398874989479, 0.61803398874989479, 4.6180339887498949 ], [ 4.6180339887498949, 0.61803398874989479, 0.61803398874989479 ], [ 0.61803398874989479, 4.6180339887498949, 0.61803398874989479 ], [ -0.61803398874989479, 0.61803398874989479, 4.6180339887498949 ], [ 4.6180339887498949, -0.61803398874989479, 0.61803398874989479 ], [ 0.61803398874989479, 4.6180339887498949, -0.61803398874989479 ], [ 0.61803398874989479, -0.61803398874989479, 4.6180339887498949 ], [ 4.6180339887498949, 0.61803398874989479, -0.61803398874989479 ], [ -0.61803398874989479, 4.6180339887498949, 0.61803398874989479 ], [ 0.61803398874989479, 0.61803398874989479, -4.6180339887498949 ], [ -4.6180339887498949, 0.61803398874989479, 0.61803398874989479 ], [ 0.61803398874989479, -4.6180339887498949, 0.61803398874989479 ], [ -0.61803398874989479, -0.61803398874989479, 4.6180339887498949 ], [ 4.6180339887498949, -0.61803398874989479, -0.61803398874989479 ], [ -0.61803398874989479, 4.6180339887498949, -0.61803398874989479 ], [ -0.61803398874989479, 0.61803398874989479, -4.6180339887498949 ], [ -4.6180339887498949, -0.61803398874989479, 0.61803398874989479 ], [ 0.61803398874989479, -4.6180339887498949, -0.61803398874989479 ], [ 0.61803398874989479, -0.61803398874989479, -4.6180339887498949 ], [ -4.6180339887498949, 0.61803398874989479, -0.61803398874989479 ], [ -0.61803398874989479, -4.6180339887498949, 0.61803398874989479 ], [ -0.61803398874989479, -0.61803398874989479, -4.6180339887498949 ], [ -4.6180339887498949, -0.61803398874989479, -0.61803398874989479 ], [ -0.61803398874989479, -4.6180339887498949, -0.61803398874989479 ], [ 1.2360679774997896, 1.6180339887498949, 4.2360679774997898 ], [ 4.2360679774997898, 1.2360679774997896, 1.6180339887498949 ], [ 1.6180339887498949, 4.2360679774997898, 1.2360679774997896 ], [ -1.2360679774997896, 1.6180339887498949, 4.2360679774997898 ], [ 4.2360679774997898, -1.2360679774997896, 1.6180339887498949 ], [ 1.6180339887498949, 4.2360679774997898, -1.2360679774997896 ], [ 1.2360679774997896, -1.6180339887498949, 4.2360679774997898 ], [ 4.2360679774997898, 1.2360679774997896, -1.6180339887498949 ], [ -1.6180339887498949, 4.2360679774997898, 1.2360679774997896 ], [ 1.2360679774997896, 1.6180339887498949, -4.2360679774997898 ], [ -4.2360679774997898, 1.2360679774997896, 1.6180339887498949 ], [ 1.6180339887498949, -4.2360679774997898, 1.2360679774997896 ], [ -1.2360679774997896, -1.6180339887498949, 4.2360679774997898 ], [ 4.2360679774997898, -1.2360679774997896, -1.6180339887498949 ], [ -1.6180339887498949, 4.2360679774997898, -1.2360679774997896 ], [ -1.2360679774997896, 1.6180339887498949, -4.2360679774997898 ], [ -4.2360679774997898, -1.2360679774997896, 1.6180339887498949 ], [ 1.6180339887498949, -4.2360679774997898, -1.2360679774997896 ], [ 1.2360679774997896, -1.6180339887498949, -4.2360679774997898 ], [ -4.2360679774997898, 1.2360679774997896, -1.6180339887498949 ], [ -1.6180339887498949, -4.2360679774997898, 1.2360679774997896 ], [ -1.2360679774997896, -1.6180339887498949, -4.2360679774997898 ], [ -4.2360679774997898, -1.2360679774997896, -1.6180339887498949 ], [ -1.6180339887498949, -4.2360679774997898, -1.2360679774997896 ], [ 0.61803398874989479, 2.6180339887498949, 3.8541019662496847 ], [ 3.8541019662496847, 0.61803398874989479, 2.6180339887498949 ], [ 2.6180339887498949, 3.8541019662496847, 0.61803398874989479 ], [ -0.61803398874989479, 2.6180339887498949, 3.8541019662496847 ], [ 3.8541019662496847, -0.61803398874989479, 2.6180339887498949 ], [ 2.6180339887498949, 3.8541019662496847, -0.61803398874989479 ], [ 0.61803398874989479, -2.6180339887498949, 3.8541019662496847 ], [ 3.8541019662496847, 0.61803398874989479, -2.6180339887498949 ], [ -2.6180339887498949, 3.8541019662496847, 0.61803398874989479 ], [ 0.61803398874989479, 2.6180339887498949, -3.8541019662496847 ], [ -3.8541019662496847, 0.61803398874989479, 2.6180339887498949 ], [ 2.6180339887498949, -3.8541019662496847, 0.61803398874989479 ], [ -0.61803398874989479, -2.6180339887498949, 3.8541019662496847 ], [ 3.8541019662496847, -0.61803398874989479, -2.6180339887498949 ], [ -2.6180339887498949, 3.8541019662496847, -0.61803398874989479 ], [ -0.61803398874989479, 2.6180339887498949, -3.8541019662496847 ], [ -3.8541019662496847, -0.61803398874989479, 2.6180339887498949 ], [ 2.6180339887498949, -3.8541019662496847, -0.61803398874989479 ], [ 0.61803398874989479, -2.6180339887498949, -3.8541019662496847 ], [ -3.8541019662496847, 0.61803398874989479, -2.6180339887498949 ], [ -2.6180339887498949, -3.8541019662496847, 0.61803398874989479 ], [ -0.61803398874989479, -2.6180339887498949, -3.8541019662496847 ], [ -3.8541019662496847, -0.61803398874989479, -2.6180339887498949 ], [ -2.6180339887498949, -3.8541019662496847, -0.61803398874989479 ], [ 2.2360679774997898, 2, 3.6180339887498949 ], [ 3.6180339887498949, 2.2360679774997898, 2 ], [ 2, 3.6180339887498949, 2.2360679774997898 ], [ -2.2360679774997898, 2, 3.6180339887498949 ], [ 3.6180339887498949, -2.2360679774997898, 2 ], [ 2, 3.6180339887498949, -2.2360679774997898 ], [ 2.2360679774997898, -2, 3.6180339887498949 ], [ 3.6180339887498949, 2.2360679774997898, -2 ], [ -2, 3.6180339887498949, 2.2360679774997898 ], [ 2.2360679774997898, 2, -3.6180339887498949 ], [ -3.6180339887498949, 2.2360679774997898, 2 ], [ 2, -3.6180339887498949, 2.2360679774997898 ], [ -2.2360679774997898, -2, 3.6180339887498949 ], [ 3.6180339887498949, -2.2360679774997898, -2 ], [ -2, 3.6180339887498949, -2.2360679774997898 ], [ -2.2360679774997898, 2, -3.6180339887498949 ], [ -3.6180339887498949, -2.2360679774997898, 2 ], [ 2, -3.6180339887498949, -2.2360679774997898 ], [ 2.2360679774997898, -2, -3.6180339887498949 ], [ -3.6180339887498949, 2.2360679774997898, -2 ], [ -2, -3.6180339887498949, 2.2360679774997898 ], [ -2.2360679774997898, -2, -3.6180339887498949 ], [ -3.6180339887498949, -2.2360679774997898, -2 ], [ -2, -3.6180339887498949, -2.2360679774997898 ], [ 1.6180339887498949, 3, 3.2360679774997898 ], [ 3.2360679774997898, 1.6180339887498949, 3 ], [ 3, 3.2360679774997898, 1.6180339887498949 ], [ -1.6180339887498949, 3, 3.2360679774997898 ], [ 3.2360679774997898, -1.6180339887498949, 3 ], [ 3, 3.2360679774997898, -1.6180339887498949 ], [ 1.6180339887498949, -3, 3.2360679774997898 ], [ 3.2360679774997898, 1.6180339887498949, -3 ], [ -3, 3.2360679774997898, 1.6180339887498949 ], [ 1.6180339887498949, 3, -3.2360679774997898 ], [ -3.2360679774997898, 1.6180339887498949, 3 ], [ 3, -3.2360679774997898, 1.6180339887498949 ], [ -1.6180339887498949, -3, 3.2360679774997898 ], [ 3.2360679774997898, -1.6180339887498949, -3 ], [ -3, 3.2360679774997898, -1.6180339887498949 ], [ -1.6180339887498949, 3, -3.2360679774997898 ], [ -3.2360679774997898, -1.6180339887498949, 3 ], [ 3, -3.2360679774997898, -1.6180339887498949 ], [ 1.6180339887498949, -3, -3.2360679774997898 ], [ -3.2360679774997898, 1.6180339887498949, -3 ], [ -3, -3.2360679774997898, 1.6180339887498949 ], [ -1.6180339887498949, -3, -3.2360679774997898 ], [ -3.2360679774997898, -1.6180339887498949, -3 ], [ -3, -3.2360679774997898, -1.6180339887498949 ] ],[ [ 1, 4 ], [ 1, 7 ], [ 1, 25 ], [ 2, 5 ], [ 2, 8 ], [ 2, 26 ], [ 3, 6 ], [ 3, 9 ], [ 3, 27 ], [ 4, 13 ], [ 4, 28 ], [ 5, 14 ], [ 5, 29 ], [ 6, 15 ], [ 6, 30 ], [ 7, 13 ], [ 7, 31 ], [ 8, 14 ], [ 8, 32 ], [ 9, 15 ], [ 9, 33 ], [ 10, 16 ], [ 10, 19 ], [ 10, 34 ], [ 11, 17 ], [ 11, 20 ], [ 11, 35 ], [ 12, 18 ], [ 12, 21 ], [ 12, 36 ], [ 13, 37 ], [ 14, 38 ], [ 15, 39 ], [ 16, 22 ], [ 16, 40 ], [ 17, 23 ], [ 17, 41 ], [ 18, 24 ], [ 18, 42 ], [ 19, 22 ], [ 19, 43 ], [ 20, 23 ], [ 20, 44 ], [ 21, 24 ], [ 21, 45 ], [ 22, 46 ], [ 23, 47 ], [ 24, 48 ], [ 25, 49 ], [ 25, 73 ], [ 26, 50 ], [ 26, 74 ], [ 27, 51 ], [ 27, 75 ], [ 28, 52 ], [ 28, 76 ], [ 29, 53 ], [ 29, 77 ], [ 30, 54 ], [ 30, 78 ], [ 31, 55 ], [ 31, 79 ], [ 32, 56 ], [ 32, 80 ], [ 33, 57 ], [ 33, 81 ], [ 34, 58 ], [ 34, 82 ], [ 35, 59 ], [ 35, 83 ], [ 36, 60 ], [ 36, 84 ], [ 37, 61 ], [ 37, 85 ], [ 38, 62 ], [ 38, 86 ], [ 39, 63 ], [ 39, 87 ], [ 40, 64 ], [ 40, 88 ], [ 41, 65 ], [ 41, 89 ], [ 42, 66 ], [ 42, 90 ], [ 43, 67 ], [ 43, 91 ], [ 44, 68 ], [ 44, 92 ], [ 45, 69 ], [ 45, 93 ], [ 46, 70 ], [ 46, 94 ], [ 47, 71 ], [ 47, 95 ], [ 48, 72 ], [ 48, 96 ], [ 49, 52 ], [ 49, 97 ], [ 50, 53 ], [ 50, 98 ], [ 51, 54 ], [ 51, 99 ], [ 52, 100 ], [ 53, 101 ], [ 54, 102 ], [ 55, 61 ], [ 55, 103 ], [ 56, 62 ], [ 56, 104 ], [ 57, 63 ], [ 57, 105 ], [ 58, 64 ], [ 58, 106 ], [ 59, 65 ], [ 59, 107 ], [ 60, 66 ], [ 60, 108 ], [ 61, 109 ], [ 62, 110 ], [ 63, 111 ], [ 64, 112 ], [ 65, 113 ], [ 66, 114 ], [ 67, 70 ], [ 67, 115 ], [ 68, 71 ], [ 68, 116 ], [ 69, 72 ], [ 69, 117 ], [ 70, 118 ], [ 71, 119 ], [ 72, 120 ], [ 73, 97 ], [ 73, 98 ], [ 74, 98 ], [ 74, 99 ], [ 75, 97 ], [ 75, 99 ], [ 76, 100 ], [ 76, 107 ], [ 77, 101 ], [ 77, 108 ], [ 78, 102 ], [ 78, 106 ], [ 79, 101 ], [ 79, 103 ], [ 80, 102 ], [ 80, 104 ], [ 81, 100 ], [ 81, 105 ], [ 82, 104 ], [ 82, 106 ], [ 83, 105 ], [ 83, 107 ], [ 84, 103 ], [ 84, 108 ], [ 85, 109 ], [ 85, 113 ], [ 86, 110 ], [ 86, 114 ], [ 87, 111 ], [ 87, 112 ], [ 88, 112 ], [ 88, 116 ], [ 89, 113 ], [ 89, 117 ], [ 90, 114 ], [ 90, 115 ], [ 91, 110 ], [ 91, 115 ], [ 92, 111 ], [ 92, 116 ], [ 93, 109 ], [ 93, 117 ], [ 94, 118 ], [ 94, 119 ], [ 95, 119 ], [ 95, 120 ], [ 96, 118 ], [ 96, 120 ] ],[ [ 1, 4, 7, 13 ], [ 2, 5, 8, 14 ], [ 3, 6, 9, 15 ], [ 10, 16, 19, 22 ], [ 11, 17, 20, 23 ], [ 12, 18, 21, 24 ], [ 25, 49, 73, 97 ], [ 26, 50, 74, 98 ], [ 27, 51, 75, 99 ], [ 28, 52, 76, 100 ], [ 29, 53, 77, 101 ], [ 30, 54, 78, 102 ], [ 31, 55, 79, 103 ], [ 32, 56, 80, 104 ], [ 33, 57, 81, 105 ], [ 34, 58, 82, 106 ], [ 35, 59, 83, 107 ], [ 36, 60, 84, 108 ], [ 37, 61, 85, 109 ], [ 38, 62, 86, 110 ], [ 39, 63, 87, 111 ], [ 40, 64, 88, 112 ], [ 41, 65, 89, 113 ], [ 42, 66, 90, 114 ], [ 43, 67, 91, 115 ], [ 44, 68, 92, 116 ], [ 45, 69, 93, 117 ], [ 46, 70, 94, 118 ], [ 47, 71, 95, 119 ], [ 48, 72, 96, 120 ], [ 1, 4, 25, 28, 49, 52 ], [ 2, 5, 26, 29, 50, 53 ], [ 3, 6, 27, 30, 51, 54 ], [ 7, 13, 31, 37, 55, 61 ], [ 8, 14, 32, 38, 56, 62 ], [ 9, 15, 33, 39, 57, 63 ], [ 10, 16, 34, 40, 58, 64 ], [ 11, 17, 35, 41, 59, 65 ], [ 12, 18, 36, 42, 60, 66 ], [ 19, 22, 43, 46, 67, 70 ], [ 20, 23, 44, 47, 68, 71 ], [ 21, 24, 45, 48, 69, 72 ], [ 73, 74, 75, 97, 98, 99 ], [ 76, 81, 83, 100, 105, 107 ], [ 77, 79, 84, 101, 103, 108 ], [ 78, 80, 82, 102, 104, 106 ], [ 85, 89, 93, 109, 113, 117 ], [ 86, 90, 91, 110, 114, 115 ], [ 87, 88, 92, 111, 112, 116 ], [ 94, 95, 96, 118, 119, 120 ], [ 1, 7, 25, 31, 50, 53, 73, 79, 98, 101 ], [ 2, 8, 26, 32, 51, 54, 74, 80, 99, 102 ], [ 3, 9, 27, 33, 49, 52, 75, 81, 97, 100 ], [ 4, 13, 28, 37, 59, 65, 76, 85, 107, 113 ], [ 5, 14, 29, 38, 60, 66, 77, 86, 108, 114 ], [ 6, 15, 30, 39, 58, 64, 78, 87, 106, 112 ], [ 10, 19, 34, 43, 56, 62, 82, 91, 104, 110 ], [ 11, 20, 35, 44, 57, 63, 83, 92, 105, 111 ], [ 12, 21, 36, 45, 55, 61, 84, 93, 103, 109 ], [ 16, 22, 40, 46, 68, 71, 88, 94, 116, 119 ], [ 17, 23, 41, 47, 69, 72, 89, 95, 117, 120 ], [ 18, 24, 42, 48, 67, 70, 90, 96, 115, 118 ] ])

# ####SnubDodecahedron
#  SnubDodecahedron=Polyhedron([ [ 0.14883199999999999, 0.168576, 0.94315099999999996 ], [ -0.88953400000000005, -0.28920299999999999, -0.25533 ], [ 0.74070199999999997, -0.38118600000000002, -0.496145 ], [ -0.78534800000000005, -0.56196500000000005, 0.086754100000000001 ], [ -0.086754100000000001, 0.78534800000000005, -0.56196500000000005 ], [ 0.56196500000000005, 0.086754100000000001, 0.78534800000000005 ], [ 0.496145, -0.74070199999999997, -0.38118600000000002 ], [ -0.94315099999999996, -0.14883199999999999, 0.168576 ], [ 0.25533, 0.88953400000000005, -0.28920299999999999 ], [ 0.28920299999999999, -0.25533, 0.88953400000000005 ], [ -0.168576, 0.94315099999999996, -0.14883199999999999 ], [ 0.38118600000000002, -0.496145, -0.74070199999999997 ], [ 0.14883199999999999, -0.168576, -0.94315099999999996 ], [ 0.74070199999999997, 0.38118600000000002, 0.496145 ], [ -0.88953400000000005, 0.28920299999999999, 0.25533 ], [ -0.78534800000000005, 0.56196500000000005, -0.086754100000000001 ], [ 0.56196500000000005, -0.086754100000000001, -0.78534800000000005 ], [ -0.086754100000000001, -0.78534800000000005, 0.56196500000000005 ], [ 0.28920299999999999, 0.25533, -0.88953400000000005 ], [ -0.168576, -0.94315099999999996, 0.14883199999999999 ], [ 0.38118600000000002, 0.496145, 0.74070199999999997 ], [ 0.496145, 0.74070199999999997, 0.38118600000000002 ], [ -0.94315099999999996, 0.14883199999999999, -0.168576 ], [ 0.25533, -0.88953400000000005, 0.28920299999999999 ], [ 0.63651599999999997, 0.65394799999999997, -0.327569 ], [ 0.327569, -0.63651599999999997, 0.65394799999999997 ], [ -0.65394799999999997, -0.327569, -0.63651599999999997 ], [ 0.63651599999999997, -0.65394799999999997, 0.327569 ], [ -0.65394799999999997, 0.327569, 0.63651599999999997 ], [ 0.327569, 0.63651599999999997, -0.65394799999999997 ], [ -0.327569, -0.63651599999999997, -0.65394799999999997 ], [ 0.65394799999999997, -0.327569, 0.63651599999999997 ], [ -0.63651599999999997, 0.65394799999999997, 0.327569 ], [ 0.65394799999999997, 0.327569, -0.63651599999999997 ], [ -0.327569, 0.63651599999999997, 0.65394799999999997 ], [ -0.63651599999999997, -0.65394799999999997, -0.327569 ], [ -0.25533, 0.88953400000000005, 0.28920299999999999 ], [ 0.94315099999999996, -0.14883199999999999, -0.168576 ], [ -0.496145, -0.74070199999999997, 0.38118600000000002 ], [ -0.38118600000000002, -0.496145, 0.74070199999999997 ], [ 0.168576, 0.94315099999999996, 0.14883199999999999 ], [ -0.28920299999999999, -0.25533, -0.88953400000000005 ], [ 0.086754100000000001, 0.78534800000000005, 0.56196500000000005 ], [ -0.56196500000000005, 0.086754100000000001, -0.78534800000000005 ], [ 0.78534800000000005, -0.56196500000000005, -0.086754100000000001 ], [ 0.88953400000000005, -0.28920299999999999, 0.25533 ], [ -0.74070199999999997, -0.38118600000000002, 0.496145 ], [ -0.14883199999999999, 0.168576, -0.94315099999999996 ], [ -0.38118600000000002, 0.496145, -0.74070199999999997 ], [ 0.168576, -0.94315099999999996, -0.14883199999999999 ], [ -0.28920299999999999, 0.25533, 0.88953400000000005 ], [ -0.25533, -0.88953400000000005, -0.28920299999999999 ], [ 0.94315099999999996, 0.14883199999999999, 0.168576 ], [ -0.496145, 0.74070199999999997, -0.38118600000000002 ], [ -0.56196500000000005, -0.086754100000000001, 0.78534800000000005 ], [ 0.086754100000000001, -0.78534800000000005, -0.56196500000000005 ], [ 0.78534800000000005, 0.56196500000000005, 0.086754100000000001 ], [ -0.74070199999999997, 0.38118600000000002, -0.496145 ], [ 0.88953400000000005, 0.28920299999999999, -0.25533 ], [ -0.14883199999999999, -0.168576, 0.94315099999999996 ] ],[ [ 1, 6 ], [ 1, 10 ], [ 1, 21 ], [ 1, 51 ], [ 1, 60 ], [ 2, 4 ], [ 2, 8 ], [ 2, 23 ], [ 2, 27 ], [ 2, 36 ], [ 3, 7 ], [ 3, 12 ], [ 3, 17 ], [ 3, 38 ], [ 3, 45 ], [ 4, 8 ], [ 4, 36 ], [ 4, 39 ], [ 4, 47 ], [ 5, 9 ], [ 5, 11 ], [ 5, 30 ], [ 5, 49 ], [ 5, 54 ], [ 6, 10 ], [ 6, 14 ], [ 6, 21 ], [ 6, 32 ], [ 7, 12 ], [ 7, 45 ], [ 7, 50 ], [ 7, 56 ], [ 8, 15 ], [ 8, 23 ], [ 8, 47 ], [ 9, 11 ], [ 9, 25 ], [ 9, 30 ], [ 9, 41 ], [ 10, 26 ], [ 10, 32 ], [ 10, 60 ], [ 11, 37 ], [ 11, 41 ], [ 11, 54 ], [ 12, 13 ], [ 12, 17 ], [ 12, 56 ], [ 13, 17 ], [ 13, 19 ], [ 13, 42 ], [ 13, 48 ], [ 14, 21 ], [ 14, 22 ], [ 14, 53 ], [ 14, 57 ], [ 15, 16 ], [ 15, 23 ], [ 15, 29 ], [ 15, 33 ], [ 16, 23 ], [ 16, 33 ], [ 16, 54 ], [ 16, 58 ], [ 17, 19 ], [ 17, 34 ], [ 18, 20 ], [ 18, 24 ], [ 18, 26 ], [ 18, 39 ], [ 18, 40 ], [ 19, 30 ], [ 19, 34 ], [ 19, 48 ], [ 20, 24 ], [ 20, 39 ], [ 20, 50 ], [ 20, 52 ], [ 21, 22 ], [ 21, 43 ], [ 22, 41 ], [ 22, 43 ], [ 22, 57 ], [ 23, 58 ], [ 24, 26 ], [ 24, 28 ], [ 24, 50 ], [ 25, 30 ], [ 25, 34 ], [ 25, 57 ], [ 25, 59 ], [ 26, 28 ], [ 26, 32 ], [ 27, 31 ], [ 27, 36 ], [ 27, 42 ], [ 27, 44 ], [ 28, 32 ], [ 28, 45 ], [ 28, 46 ], [ 29, 33 ], [ 29, 35 ], [ 29, 51 ], [ 29, 55 ], [ 30, 34 ], [ 31, 36 ], [ 31, 42 ], [ 31, 52 ], [ 31, 56 ], [ 32, 46 ], [ 33, 35 ], [ 33, 37 ], [ 34, 59 ], [ 35, 37 ], [ 35, 43 ], [ 35, 51 ], [ 36, 52 ], [ 37, 41 ], [ 37, 43 ], [ 38, 45 ], [ 38, 46 ], [ 38, 53 ], [ 38, 59 ], [ 39, 40 ], [ 39, 47 ], [ 40, 47 ], [ 40, 55 ], [ 40, 60 ], [ 41, 43 ], [ 42, 44 ], [ 42, 48 ], [ 44, 48 ], [ 44, 49 ], [ 44, 58 ], [ 45, 46 ], [ 46, 53 ], [ 47, 55 ], [ 48, 49 ], [ 49, 54 ], [ 49, 58 ], [ 50, 52 ], [ 50, 56 ], [ 51, 55 ], [ 51, 60 ], [ 52, 56 ], [ 53, 57 ], [ 53, 59 ], [ 54, 58 ], [ 55, 60 ], [ 57, 59 ] ],[ [ 53, 57, 59 ], [ 51, 55, 60 ], [ 50, 52, 56 ], [ 49, 54, 58 ], [ 44, 48, 49 ], [ 44, 49, 58 ], [ 42, 44, 48 ], [ 40, 47, 55 ], [ 40, 55, 60 ], [ 39, 40, 47 ], [ 38, 45, 46 ], [ 38, 46, 53 ], [ 38, 53, 59 ], [ 37, 41, 43 ], [ 35, 37, 43 ], [ 33, 35, 37 ], [ 31, 36, 52 ], [ 31, 52, 56 ], [ 29, 33, 35 ], [ 29, 35, 51 ], [ 29, 51, 55 ], [ 28, 32, 46 ], [ 28, 45, 46 ], [ 27, 31, 36 ], [ 27, 31, 42 ], [ 27, 42, 44 ], [ 26, 28, 32 ], [ 25, 30, 34 ], [ 25, 34, 59 ], [ 25, 57, 59 ], [ 24, 26, 28 ], [ 22, 41, 43 ], [ 21, 22, 43 ], [ 20, 24, 50 ], [ 20, 50, 52 ], [ 19, 30, 34 ], [ 18, 20, 24 ], [ 18, 20, 39 ], [ 18, 24, 26 ], [ 18, 39, 40 ], [ 17, 19, 34 ], [ 16, 23, 58 ], [ 16, 54, 58 ], [ 15, 16, 23 ], [ 15, 16, 33 ], [ 15, 29, 33 ], [ 14, 21, 22 ], [ 14, 22, 57 ], [ 14, 53, 57 ], [ 13, 17, 19 ], [ 13, 19, 48 ], [ 13, 42, 48 ], [ 12, 13, 17 ], [ 12, 13, 31, 42, 56 ], [ 11, 37, 41 ], [ 11, 16, 33, 37, 54 ], [ 10, 26, 32 ], [ 10, 18, 26, 40, 60 ], [ 9, 25, 30 ], [ 9, 11, 41 ], [ 9, 22, 25, 41, 57 ], [ 8, 15, 23 ], [ 8, 15, 29, 47, 55 ], [ 7, 12, 56 ], [ 7, 24, 28, 45, 50 ], [ 7, 50, 56 ], [ 6, 14, 21 ], [ 6, 10, 32 ], [ 6, 14, 32, 46, 53 ], [ 5, 9, 11 ], [ 5, 11, 54 ], [ 5, 9, 30 ], [ 5, 19, 30, 48, 49 ], [ 5, 49, 54 ], [ 4, 8, 47 ], [ 4, 20, 36, 39, 52 ], [ 4, 39, 47 ], [ 3, 7, 12 ], [ 3, 7, 45 ], [ 3, 12, 17 ], [ 3, 17, 34, 38, 59 ], [ 3, 38, 45 ], [ 2, 4, 8 ], [ 2, 4, 36 ], [ 2, 8, 23 ], [ 2, 23, 27, 44, 58 ], [ 2, 27, 36 ], [ 1, 6, 10 ], [ 1, 10, 60 ], [ 1, 6, 21 ], [ 1, 21, 35, 43, 51 ], [ 1, 51, 60 ] ]; atol = 1e-6)

# #### TriangularPrism
#  TriangularPrism=Polyhedron([ [ 0, 0, 0 ], [ 1, 0, 0 ], [ 1/2, 0.8660254037844386, 0 ], [ 1/2, 0.8660254037844386, 1 ], [ 0, 0, 1 ], [ 1, 0, 1 ] ],[ [ 1, 2 ], [ 1, 3 ], [ 1, 5 ], [ 2, 3 ], [ 2, 6 ], [ 3, 4 ], [ 4, 5 ], [ 4, 6 ], [ 5, 6 ] ],[ [ 1, 2, 3 ], [ 1, 3, 4, 5 ], [ 2, 3, 4, 6 ], [ 1, 2, 5, 6 ], [ 4, 5, 6 ] ])

# #### HexagonalPrism
#  HexagonalPrism=Polyhedron([ [ 0.5, 0.8660254037844386, 0. ], [ -0.5, 0.8660254037844386, 0. ], [ -1., 0., 0. ], [ -0.5, -0.8660254037844386, 0. ], [ 1/2, -0.8660254037844386, 0 ], [ 1., 0, 0 ], [ 1/2, 0.8660254037844386, 1 ], [ -1/2, 0.8660254037844386, 1 ], [ -1., 0., 1 ], [ -1/2, -0.8660254037844386, 1 ], [ 1/2, -0.8660254037844386, 1 ], [ 1., 0, 1 ] ],[ [ 1, 2 ], [ 1, 6 ], [ 1, 7 ], [ 2, 3 ], [ 2, 8 ], [ 3, 4 ], [ 3, 9 ], [ 4, 5 ], [ 4, 10 ], [ 5, 6 ], [ 5, 11 ], [ 6, 12 ], [ 7, 8 ], [ 7, 12 ], [ 8, 9 ], [ 9, 10 ], [ 10, 11 ], [ 11, 12 ] ],[ [ 1, 2, 3, 4, 5, 6 ], [ 7, 8, 9, 10, 11, 12 ], [ 1, 2, 7, 8 ], [ 2, 3, 8, 9 ], [ 3, 4, 9, 10 ], [ 4, 5, 10, 11 ], [ 5, 6, 11, 12 ], [ 1, 6, 7, 12 ] ])

# ####Johnson
#  Johnson=Polyhedron([ [ 1/2, 0, -0.8660254037844386 ], [ 1/2, 1, -0.8660254037844386 ], [ 0, 1, 0 ], [ 0, 0, 0 ], [ 1, 1, 0 ], [ 1, 0, 0 ], [ 0, 1/2, 0.8660254037844386 ], [ 1, 1/2, 0.8660254037844386 ] ],[ [ 1, 2 ], [ 1, 4 ], [ 1, 6 ], [ 2, 3 ], [ 2, 5 ], [ 3, 4 ], [ 3, 5 ], [ 3, 7 ], [ 4, 6 ], [ 4, 7 ], [ 5, 6 ], [ 5, 8 ], [ 6, 8 ], [ 7, 8 ] ],[ [ 1, 2, 3, 4 ], [ 2, 3, 5 ], [ 1, 2, 5, 6 ], [ 1, 4, 6 ], [ 3, 4, 7 ], [ 3, 5, 7, 8 ], [ 5, 6, 8 ], [ 4, 6, 7, 8 ] ])

# #### EloncatedDodecahedron
#  EloncatedDodecahedron=Polyhedron([ [ 0.8660254037844386, 0, 3/2 ], [ 1.7320508075688772, 0, 1 ], [ 1.7320508075688772, 0, 0 ], [ 0.8660254037844386, 0, -1/2 ], [ 0, 0, 0 ], [ 0, 0, 1 ], [ 1.7320508075688772, 0.8660254037844386, 3/2 ], [ 1.7320508075688772, 1.7320508075688772, 1 ], [ 1.7320508075688772, 1.7320508075688772, 0 ], [ 1.7320508075688772, 0.8660254037844386, -1/2 ], [ 0.8660254037844386, 1.7320508075688772, 3/2 ], [ 0, 1.7320508075688772, 1 ], [ 0, 1.7320508075688772, 0 ], [ 0.8660254037844386, 1.7320508075688772, -1/2 ], [ 0, 0.8660254037844386, 3/2 ], [ 0, 0.8660254037844386, -1/2 ], [ 0.8660254037844386, 0.8660254037844386, -1 ], [ 0.8660254037844386, 0.8660254037844386, 2 ] ],[ [ 1, 2 ], [ 1, 6 ], [ 1, 18 ], [ 2, 3 ], [ 2, 7 ], [ 3, 4 ], [ 3, 10 ], [ 4, 5 ], [ 4, 17 ], [ 5, 6 ], [ 5, 16 ], [ 6, 15 ], [ 7, 8 ], [ 7, 18 ], [ 8, 9 ], [ 8, 11 ], [ 9, 10 ], [ 9, 14 ], [ 10, 17 ], [ 11, 12 ], [ 11, 18 ], [ 12, 13 ], [ 12, 15 ], [ 13, 14 ], [ 13, 16 ], [ 14, 17 ], [ 15, 18 ], [ 16, 17 ] ],[ [ 1, 2, 3, 4, 5, 6 ], [ 2, 3, 7, 8, 9, 10 ], [ 8, 9, 11, 12, 13, 14 ], [ 5, 6, 12, 13, 15, 16 ], [ 1, 6, 15, 18 ], [ 1, 2, 7, 18 ], [ 7, 8, 11, 18 ], [ 11, 12, 15, 18 ], [ 4, 5, 16, 17 ], [ 3, 4, 10, 17 ], [ 9, 10, 14, 17 ], [ 13, 14, 16, 17 ] ])


##################################################################################
################# Bricard Octahedron 
##################################################################################
"""
    bricard_octahedron(n::Integer)

Return the bricard octahedron. n determines, what symmetry is used to calculate the position of the second pole. 1, 2 = plane reflection, 3 = line reflection
"""
function bricard_octahedron(n::Integer)
  equator = [[-1, 0, 0], [0, -1, -1], [1, 0, 0], [0, 1, -1]]
  north = [-0.5, 0.5, 2]

  if n == 1
    aff = rigidmap(hcat(equator..., [0, 0, 1]), hcat(equator[3], equator[4], equator[1], equator[2], [0, 0, 1]))
  elseif n == 2
    aff = rigidmap(hcat(equator..., [0, 0, 1]), hcat(equator[1], equator[4], equator[3], equator[2], [0, 0, 1]))
  elseif n == 3
    aff = rigidmap(hcat(equator..., [0, 0, 1]), hcat(equator[3], equator[2], equator[1], equator[4], [0, 0, 1]))
  else
    @error "n needs to be either 1,2 or 3, but it is $n."
  end

  south = aff(north)

  verts = vcat(equator, [north, south])
  edges = [[1, 2], [2, 3], [3, 4], [4, 1], [1, 5], [2, 5], [3, 5], [4, 5], [1, 6], [2, 6], [3, 6], [4, 6]]
  faces = [[1, 2, 5], [2, 3, 5], [3, 4, 5], [1, 4, 5], [1, 2, 6], [2, 3, 6], [3, 4, 6], [1, 4, 6]]

  return Polyhedron(verts, edges, faces)
end

# ##################################################################################
# ################# further example
# ##################################################################################

# #example for ramified complex

# function star(n::Int)

#     vertices=[[0.0,0.0,0.5],[0.0,0.0,-0.5]]

#     facets=[]
#     edges=[[1,2]]
#     for i in 1:n
# 	push!(edges,[1,i+2])
# 	push!(edges,[2,i+2])
# 	push!(facets,[1,2,i+2])
#       push!(vertices,[sqrt(3.0)/2.0*cos(2*3.1415/n),sqrt(3.0)/2*cos(2.0*3.1415/n),0.0])
#     end
#    return Polyhedron(vertices,edges,facets) 
# end


# function Prisma(n::Int)
#   #vertices 
#   t=1/sqrt((cos(3.14159*2*1/n)-cos(3.14159*2*2/n))^2+
#            (sin(3.14159*2*1/n)-sin(3.14159*2*2/n))^2);
#   vertices=Vector{Vector{Float64}}(undef,0)
#   vertices=map(i->[t*cos(3.14159*2*i/n),t*sin(3.14159*2*i/n),0.0],1:n)
#   for i in 1:n 
#     push!(vertices,[t*cos(3.14159*2*i/n),t*sin(3.14159*2*i/n),1.0])
#   end
#   #facets
#   facets=Vector{Vector{Int}}(undef,0)
#   for i in 1:n-1
#     push!(facets,[i,i+1,n+i+1,n+i])
#   end 
#   push!(facets,[n,1,n+1,2*n]) 
#   push!(facets,Vector((1:n)))
#   push!(facets,Vector((n+1:2*n)))
#   #edges
#   edges=Vector{Vector{Int}}(undef,0)
#   for i in 1:n-1
#     push!(edges,[i,i+1])
#     push!(edges,[n+i,n+i+1])
#     push!(edges,[i,n+i])
#   end
#   push!(edges,[1,n])
#   push!(edges,[n+1,2*n])
#   push!(edges,[n,2*n])
#   return Polyhedron(vertices,edges,facets) 
# end;


