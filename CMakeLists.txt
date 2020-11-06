# 
# INCLUDE_DIRECTORIES(EquationSystems)
# INCLUDE_DIRECTORIES(AdvectionTerms)

SET(WeakMHDSolverSource 
#         ./EquationSystems/CoupledLinearNS.cpp
#         ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
	./EquationSystems/IncNavierStokes.cpp
        ./EquationSystems/VelocityCorrectionScheme.cpp
#         ./EquationSystems/VCSMapping.cpp
	./EquationSystems/Extrapolate.cpp
        ./EquationSystems/StandardExtrapolate.cpp
#        ./EquationSystems/MappingExtrapolate.cpp
#        ./EquationSystems/SubSteppingExtrapolate.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
        ./AdvectionTerms/LinearisedAdvection.cpp
        ./AdvectionTerms/NavierStokesAdvection.cpp
        ./AdvectionTerms/NonlinearDisturbanceAdvection.cpp
#        ./AdvectionTerms/SkewSymmetricAdvection.cpp
#        ./AdvectionTerms/AlternateSkewAdvection.cpp
#        ./AdvectionTerms/NoAdvection.cpp
#        ./Filters/FilterEnergy.cpp
#        ./Filters/FilterReynoldsStresses.cpp
#        ./Filters/FilterMovingBody.cpp
#        ./Forcing/ForcingMovingBody.cpp	
#        ./Forcing/ForcingStabilityCoupledLNS.cpp
       ./SolverUtils/Driver.cpp
       ./SolverUtils/WMHDDriverArpack.cpp
       ./SolverUtils/WMHDDriverArnoldi.cpp
       ./SolverUtils/WMHDDriverModifiedArnoldi.cpp
       ./WeakMHDSolver.cpp 
)
    
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})



IF(USE_MPI)
    SET(WeakMHDSolverNAME "WeakMHDSolverMPI.out")
ELSE(USE_MPI)
    SET(WeakMHDSolverNAME "WeakMHDSolver.out")
ENDIF(USE_MPI)

ADD_EXECUTABLE(${WeakMHDSolverNAME} ${WeakMHDSolverSource})
TARGET_LINK_LIBRARIES(${WeakMHDSolverNAME} ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})

 
