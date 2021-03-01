    autoPtr<heterogeneousRadiationModel> radiation
    (
        heterogeneousRadiationModel::New(thermo.T(),porosityF,pyrolysisZone.surf(),pyrolysisZone.T())
    );
