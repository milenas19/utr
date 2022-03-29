/*
utr - Geant4 simulation of the UTR at HIGS
Copyright (C) 2017 the developing team (see README.md)

This file is part of utr.

utr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

utr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with utr.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
	G4VPhysicalVolume *Construct() override;
	void ConstructSDandField() override;

	void print_info() const;

private:
	G4double World_x;
	G4double World_y;
	G4double World_z;
	G4LogicalVolume *World_logical;

	void ConstructCollimator(G4ThreeVector local_coordinates);


};

