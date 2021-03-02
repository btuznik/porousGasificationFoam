if (adjustTimeStep and solidChemistryTimeStepControl)
{

    const scalar TsolidChemistry = pyrolysisZone.maxTime();
    const scalar dt0 = runTime.deltaTValue();

    if (min(TsolidChemistry,dt0) == TsolidChemistry)
    {
        Info << "timestep corrected for heterogeneous reactions to: " << TsolidChemistry << "[s]" << endl;
    }

    runTime.setDeltaT
    (
        min
        (
            TsolidChemistry,
            dt0
        )
    );

}
