/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-thermoMixture

Description

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "dictionary.H"
#include "mutation++.h"
#include "perfectGas.H"
#include "rrhoThermo.H"
#include "specie.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[]) {
  typedef rrhoThermo<perfectGas<specie>> ThermoType;

  dictionary dict(IFstream("thermoDict")());

  ThermoType t1("rrhoThermo", dict.subDict("air"));

  const scalar cp = t1.Cp(1e5, 300);

  Info << "t1.Cp(300, 1e5) = " << cp << " [J/kmol/K]" << endl;

  Info << "t1.W() =  " << t1.W() << endl;

  Info << "End" << nl << endl;

  return 0;
}

// ************************************************************************* //
