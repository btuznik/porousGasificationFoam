# TODO i pytania

## Pytania

* Co zrobić z DpDt
* biomassGasificationFoam.C: co z upadateChemistryTimeStep.H,
  co z if(nOuterCorr != 0), co z tym, że YEqn i hsEqn są poza pimple loopem
* UEqn -- nie można usunac U z urbulence->divDevTau(U)
* Czy my w ogóle mamy klase heterogenousMeanTemp?
  
## Zadania

* Dograc cale heterogenosRadiationModel -- glownie Sh i Shs (fe4.1)
* Usunąć pole p i sprawdzić czy działa bez. (fe4.1)
* ogarnąć co w ogóle z plikie setMultiRegionDeltaT.H
* Info w YEqn.H ogarnąć (fe4.1)
* Zdecydować czy dodać klase heterogenousMeanTemp (fe4.1)
* Co z dictionary coeffs w heterogenosRadiationModel.C -- czy zostawić tak jak
  czy jako null? (fe4.1)
* Sprawdzic czy cala radiacja jest odkomentowana

* ConstSolidThermo i ExponentialSolidThermo są nadal nie zmienione (fe4.1)
* Chemia -- zmienić (fe.41)

* Przelecieć po wszystkich TODO

## Co ostatnio skończyłem

Porównałem wersję cleanVersionDist  z RecentVersion. Na razie doszedłem do
volPyrolysis (udało się ją skończyć i się kompiluje). biomass też się kompiluje.
