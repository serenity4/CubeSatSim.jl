#= Component types =#

abstract type PowerConsumer end
abstract type PowerGenerator end

abstract type Actuator <: PowerConsumer end
abstract type Sensor <: PowerConsumer end

abstract type Antenna <: PowerConsumer end
abstract type AttitudeActuator <: Actuator end
abstract type AttitudeSensor <: Sensor end
abstract type OnBoardComputer <: PowerConsumer end
abstract type Payload <: PowerConsumer end

moments(::AttitudeActuator) = error("Not implemented by user")
power_generated(::PowerGenerator) = error("Not implemented by user")
power_consumption(::PowerConsumer) = error("Not implemented by user")
frequency(::Antenna) = error("Not implemented by user")
bandwidth(::Antenna) = error("Not implemented by user")

#= Perturbation types =#

abstract type PerturbationModel

"Aerodynamic model, containing forces and moments."
abstract type AerodynamicModel <: PerturbationModel end
abstract type MagneticModel <: PerturbationModel end
abstract type GravityModel <: PerturbationModel end

forces(::AerodynamicModel) = error("Not implemented by user")
moments(::AerodynamicModel) = error("Not implemented by user")
