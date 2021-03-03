#ifndef READCHEMISTRYTIMECONTROLS_H
#define READCHEMISTRYTIMECONTROLS_H

Switch gasChemistryTimeStepControl
(
    mesh.lookupObject<dictionary>
        ("chemistryProperties").lookupOrDefault("gasChemistryTimeStepControl", false)
);

Switch solidChemistryTimeStepControl
(
    mesh.lookupObject<dictionary>
        ("chemistryProperties").lookupOrDefault("solidChemistryTimeStepControl", false)
);

#endif // READCHEMISTRYTIMECONTROLS_H
