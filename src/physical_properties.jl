using UnPack

abstract type AbstractVolumeExpansion
end 

#------------ Single Component Structs ---------------

mutable struct CpPoly{T}
    A::T
    B::T
    C::T
end

mutable struct CpMixingRule{T <: Real}
    Cp1::T
    Cp2::T
    Cp3::T
    Cp4::T  
end

mutable struct HallvardVolumeExpansionPoly{T} <: AbstractVolumeExpansion
    a1::T 
    a2::T
    a3::T
end

mutable struct VᴱPoly{T} 
    k1::T
    k2::T
    k3::T 
    k4::T
end

mutable struct MEAViscosityPoly{T}
    b1::T
    b2::T
    b3::T
end

mutable struct H₂OViscosityPoly{T}
    A::T
    B::T
    C::T
    D::T
end

mutable struct DiffusionPoly{T}
    A::T
    B::T
    C::T
end

mutable struct PureDensityPoly{T}
    d1::T
    d2::T
    d3::T
end

mutable struct ViscosityDeviationPoly{T}
    l1::T
    l2::T
    l3::T
    l4::T
end

mutable struct ViscosityDeviationStarPoly{T}
    A::T
    B::T
    C::T
end

#------------ Data Structures for Fluid Properties ---------------

mutable struct FluidInfos{N <: AbstractString, Mw <: Real}
    name::N
    Mw::Mw
end

mutable struct SimpleMediumData{H <: Union{CpPoly, Nothing}, V <: Union{MEAViscosityPoly, H₂OViscosityPoly, Nothing}, D <: Union{PureDensityPoly, Nothing}}
    heat_capacity_data::H
    viscosity_data::V
    density_data::D
end

mutable struct SimpleFluidMedium{A <: FluidInfos, D <: SimpleMediumData}
    infos::A
    data::D
end


#------------ Polynomial Functions ---------------

function FluidInfos(;name::N, Mw::M) where {N, M}
    return FluidInfos{N, M}(name, Mw)
end

function SimpleMediumData(;heat_capacity_data::H, viscosity_data::V, density_data::D) where {H, V, D}
    return SimpleMediumData{H, V, D}(heat_capacity_data, viscosity_data, density_data)
end

function SimpleFluidMedium(;infos::I, data::D) where {I, D}
    return SimpleFluidMedium{I, D}(infos, data)
end

function CpPoly(;A, B, C)
    return CpPoly(A, B, C)
end

function VᴱPoly(;k1, k2, k3, k4)
    return VᴱPoly(k1, k2, k3, k4)
end 

function HallvardVolumeExpansionPoly(;a1, a2, a3)
    return HallvardVolumeExpansionPoly(a1, a2, a3) 
end

function DiffusionPoly(;A, B, C)
    return DiffusionPoly(A, B, C)
end

function PureDensityPoly(;d1, d2, d3)
    return PureDensityPoly(d1, d2, d3)
end 

function ViscosityDeviationPoly(;l1, l2, l3, l4)
    return ViscosityDeviationPoly(l1, l2, l3, l4)
end

function H₂OViscosityPoly(;A, B, C, D)
    return H₂OViscosityPoly(A, B, C, D)
end

function MEAViscosityPoly(;b1, b2, b3)
    return MEAViscosityPoly(b1, b2, b3)
end

function ViscosityDeviationStarPoly(;A, B, C)
    return ViscosityDeviationStarPoly(A, B, C)
end

#------------ Physical Property Functions ---------------

#------------ Density ---------------
function Vᴱ(excess_molar_volume::VᴱPoly, x_MEA::T, x_H₂O::T, t::T) where T <: Number 
    @unpack k1, k2, k3, k4 = excess_molar_volume
    Vᴱ_value = (k1 + k2 * t + k3 * x_MEA + k4 * x_MEA^2) * x_MEA * x_H₂O * 10^(-6)
    return Vᴱ_value
end

function w_CO₂_added(α::T, x_MEA::T, m::AbstractArray{T}) where T <: Number 
    αx_MEA = α * x_MEA
    m_MEA, m_H₂O, m_CO₂ = m
    num = αx_MEA * m_CO₂
    den = x_MEA * m_MEA + (1 - x_MEA - αx_MEA) * m_H₂O + αx_MEA * m_CO₂
    return num/den
end 

# ρ_unloaded expects mole fractions xᵢ, Mᵢ in kg/mol, ρᵢ in kg/m³, Vᴱ in m³/mol
function ρ_unloaded(Vᴱ::T, xᵢ::AbstractVector{T}, Mᵢ::AbstractVector{T}, ρᵢ::AbstractVector{T}) where T <: Number 
    xᵢMᵢ = xᵢ .* Mᵢ
    num = sum(xᵢMᵢ) # kg/mol total
    den = Vᴱ + sum(xᵢMᵢ ./ ρᵢ) # m³/mol
    ρ_value = num / den # kg/m³
    return ρ_value
end

function ρ_loaded(ρ_unloaded::T, w_CO2::T, volume_expansion::T) where T <: Number 
    num = ρ_unloaded
    den = 1 - w_CO2 * (1 - volume_expansion^3)
    return num/den 
end

function volume_expansion(α::T, x_MEA::T, volume_expansion::HallvardVolumeExpansionPoly) where T <: Number 
    @unpack a1, a2, a3 = volume_expansion
    num = a1 * x_MEA * α + a2 * x_MEA
    den  = a3 + x_MEA
    return num/den
end 

function pure_density(t::T, density_polynomial_values::PureDensityPoly) where T <: Number 
    @unpack d1, d2, d3 = density_polynomial_values
    return (d1 + d2 * t + d3 * t^2)  # kg/m³
end 

#------------ Heat Capacity ---------------
# Functions to compute the specific heat capacities sinlge components
function Cp(heat_capacity_model::CpPoly, t::T) where T <: Number 
    @unpack A, B, C =  heat_capacity_model 
    return A + B * t + C * t^2
end

# Function to compute the specific heat capacity of the liquid mixture (Cp)
function Cp_mix(mixing_rule::CpMixingRule, Cps::AbstractArray, x_MEA::T, t::T) where T <: Number 
    @unpack Cp1, Cp2, Cp3, Cp4 = mixing_rule
    (Cp_H₂O, Cp_MEA) = Cps
    return (1.0 - x_MEA) * Cp_H₂O + x_MEA * Cp_MEA + x_MEA * (1.0 - x_MEA) * (Cp1 + Cp2 * t + (Cp3 * x_MEA) / t^Cp4)
end

#------------ Diffusion ---------------
function diffusion_CO₂_MEA(diffusion_CO2_H₂O::T, diffusion_model::DiffusionPoly, μ_MEA::T, μ_H₂O::T) where T <: Number 
    @unpack A, B, C = diffusion_model
    return diffusion_CO2_H₂O * (μ_H₂O \ μ_MEA)
end

function diffusion_CO₂_H₂O(diffusion_model::DiffusionPoly, T)
    @unpack A, B, C = diffusion_model 
    return A * exp(B / T)
end

#------------ Viscosity ---------------
function η_H₂O(H₂O_viscosity_model::H₂OViscosityPoly, t::T) where T <: Number 
    @unpack A, B, C, D = H₂O_viscosity_model
    return A + B * t + C * t^2 + D * t^3
end

function η_MEA(MEA_viscosity_model::MEAViscosityPoly, t::T) where T <: Number
    @unpack b1, b2, b3 = MEA_viscosity_model
    return exp(b1 + b2 / (t + 273.15 - b3)) # Celsius to Kelvin
end

function η_deviation(viscosity_deviation_model::ViscosityDeviationPoly, t::T, x_MEA::T) where T <: Number 
    @unpack l1, l2, l3, l4 = viscosity_deviation_model
    x_H₂O = 1.0 - x_MEA
    return exp((l1 + l2 * t + l3 * t^2 + l4 * x_MEA) * x_MEA * x_H₂O)
end

function η_unloaded(viscosity_deviation::T, xᵢ::AbstractVector{T}, ηᵢ::AbstractVector{T}) where T <: Number 
    xᵢlnηᵢ = xᵢ .* log.(ηᵢ)
    return exp(log(viscosity_deviation) + sum(xᵢlnηᵢ))
end

function η_deviation_star(viscosity_deviation_star_model::ViscosityDeviationStarPoly, x_MEA::T, α::T) where T <: Number 
    @unpack A, B, C = viscosity_deviation_star_model
    num = A * x_MEA + B * α * x_MEA
    den = C + x_MEA
    return exp(num/den) / 10^3# Pa·s
end

function η_loaded(x_CO₂::T,η_star::T, η_unloaded::T) where T <: Number 
    return exp(x_CO₂ * log(η_star) + (1.0 - x_CO₂) * log(η_unloaded)) * 1e3  # Convert Pa·s to mPa·s
end

#------------ Component Information ---------------
CO₂_info = FluidInfos(name = "CO2", Mw = 44.01) # g/mol
MEA_info = FluidInfos(name = "MEA", Mw = 61.08) # g/mol
H₂O_info = FluidInfos(name = "H₂O", Mw = 18.02) # g/mol

#------------ Calculations ---------------
H₂O_heat_capacity = CpPoly(A = 4.1908, B =- 6.62e-4, C = 9.14e-6)
MEA_heat_capacity = CpPoly(A = 2.5749, B = 6.612e-3, C = - 1.9e-5)

excess_molar_volume = VᴱPoly(k1 = - 1.9210, k2 = 1.6792 * 10^(-3), k3 = - 3.0951, k4 = 3.4412)
volume_expansion_model = HallvardVolumeExpansionPoly(a1 = 0.29, a2 = 0.18, a3 = 0.66)

H₂O_density_coeffs = PureDensityPoly(d1 = 1002.3, d2 = -0.131, d3 = -0.00308) # T [°C], kg/m^3
MEA_density_coeffs = PureDensityPoly(d1 = 1023.75, d2 = -0.5575, d3 = -0.00187) # T [°C], kg/m^3

H₂O_viscosity_coeffs = H₂OViscosityPoly(A = 1.684 * 10^(-3), B = -4.264 * 10^(-5), C = 5.062 * 10^(-7), D = -2.244 * 10^(-9))
MEA_viscosity_coeffs = MEAViscosityPoly(b1 = -3.9303, b2 = 1021.8, b3 = 149.1969)
viscosity_deviation_coeffs = ViscosityDeviationPoly(l1 = 8.36, l2 = -4.664 * 10^(-2), l3 = 1.6 * 10^(-4), l4 = -4.14)
viscosity_deviation_star_coeffs = ViscosityDeviationStarPoly(A = 6.98, B = 10.48, C = 0.049)

diffusion_CO2_H₂O_coefficient = DiffusionPoly(A = 2.35 * 10^(-6), B = -2119.0, C = 0.8)

#const CO2 = SimpleFluidMedium(infos = CO2_info, data = SimpleMediumData(heat_capacity_data = CO2_heat_capacity, viscosity_data = nothing))
const MEA = SimpleFluidMedium(infos = MEA_info, data = SimpleMediumData(heat_capacity_data = MEA_heat_capacity, viscosity_data = MEA_viscosity_coeffs, density_data = MEA_density_coeffs))

# Mixtures
mixing_rule = CpMixingRule(-0.9198, 0.013969, 69.643, 1.5859)



#--------------------------------

#----- MultiComponent structs ----------


# Function to compute the density of the liquid mixture (ρ_l)
# function liquid_density(liquid_density_model::LiquidDensityPoly, densities, T)
#     @unpack A, B, C, D = liquid_density_model
#     (ρ_H2O, ρ_MEA) = densities
#     return (1.0 - C2) * ρ_H2O + C * ρ_MEA +
#            C2 * (1.0 - C2) * (A - B * T + (C * C) / T^D)
# end

# Enhancement factor (α)
# = α_expr = n_CO2 / n_MEA