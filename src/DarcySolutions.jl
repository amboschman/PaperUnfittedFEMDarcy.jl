#======================================================================# 
#                         Solution types                                #
#======================================================================#

abstract type SolutionType end

abstract type DarcySolution <: SolutionType end

#----------------------------------------------------------------------# 
#             Manufactured solutions for the Darcy Problem             # 
#----------------------------------------------------------------------#

# Case: Smooth Solenoidal Trigonometric Darcy2D
SmoothSolenoidalDarcy2D_p(x) = sin(π*x[1])-sin(π*x[2])
SmoothSolenoidalDarcy2D_f(x) = ∇(SmoothSolenoidalDarcy2D_p)(x) + VectorValue(x[1]+sin(π*x[2]),-x[2]+sin(π*x[1]))
SmoothSolenoidalDarcy2D_u(x) = SmoothSolenoidalDarcy2D_f(x) - ∇(SmoothSolenoidalDarcy2D_p)(x)
SmoothSolenoidalDarcy2D_g(x) = 0.0
SmoothSolenoidalDarcy2D_kinv = TensorValue(1.0,0.0,0.0,1.0)
SmoothSolenoidalDarcy2D_name = "SmoothSolenoidalDarcy2D"

# Case: Manufactured Pressure Robust 2-D Darcy solution (p!∈Qh,u∈Uh and divu=0)
PresRobustDarcy2D_p(x) = sin(π*x[1])-sin(π*x[2])
PresRobustDarcy2D_f(x) = ∇(PresRobustDarcy2D_p)(x) + VectorValue(x[1],-x[2])
PresRobustDarcy2D_u(x) = PresRobustDarcy2D_f(x) - ∇(PresRobustDarcy2D_p)(x)
PresRobustDarcy2D_g(x) = 0.0
PresRobustDarcy2D_kinv = TensorValue(1.0,0.0,0.0,1.0)
PresRobustDarcy2D_name = "PresRobustDarcy2D"