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

/*
 * Setup of the Sn112 NRF experiment from the 2015 DHIPS campaign
 * The purpose of this experiment was to measure the excitation strength
 * of the 2^+_1 state of 112Sn in NRF relative to the known strengths
 * of excited states in 27Al and 59Co. This was supposed to clarify the
 * experimental situation of the systematics of B(E2; 0^+_1 -> 2^+_1)
 * strengths in the tin isotopic chain.
 * 
 * An alternative would have been to determine two B(E2) strengths
 * of the tin isotopes 112Sn and 116Sn simultaneously. A sandwich
 * target which would be necessary for this alternative experiment
 * has also been implemented.
*/

#include "DetectorConstruction.hh"

// Materials
#include "G4Material.hh"
#include "G4NistManager.hh"

// Geometry
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

// Sensitive Detectors
#include "G4SDManager.hh"
#include "EnergyDepositionSD.hh"
#include "ParticleSD.hh"

// Units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


// #include "G4UnitsTable.hh"
#include "utrConfig.h"
#include <iostream>


G4VPhysicalVolume *DetectorConstruction::Construct() {

	/***************** Define Lengths ************/ 

	const double detector_to_bremstarget = 1180.*mm; // Position guesstimated
	const double bremstarget_thickness = 2.5*mm;
	const double bremstarget_edge_length = 10*mm;
	const double detector_radius = 5*mm;
	const double detector_length = 1*mm;
	const double world_buffer_length = 10*mm;

	World_x = (detector_radius + world_buffer_length) * 2;
	World_y = (detector_radius + world_buffer_length) * 2;
	World_z = detector_to_bremstarget + bremstarget_thickness + detector_length + 2 * world_buffer_length;


	/***************** Define Materials ************/

	G4NistManager *nist = G4NistManager::Instance();
	G4Material *vacuum = nist->FindOrBuildMaterial("G4_Galactic");  	//Vacuum
	G4Material *gold = nist->FindOrBuildMaterial("G4_Au");              //Bremstarget
	G4Material *detector_material = nist->FindOrBuildMaterial("G4_Pb");
	

	/******************** WORLD ******************/
	G4Box *World_solid = new G4Box("World_solid", World_x * 0.5, World_y * 0.5, World_z * 0.5);

	G4LogicalVolume *World_logical = new G4LogicalVolume(World_solid, vacuum, "World_logical", 0, 0, 0);

	//Visualisierung der Welt (Farbe)
	G4VisAttributes *world_vis = new G4VisAttributes(true, G4Color::Red());
	world_vis->SetForceWireframe(true);

	World_logical->SetVisAttributes(world_vis);

	G4VPhysicalVolume *World_physical = new G4PVPlacement(0, G4ThreeVector(), World_logical, "World", 0, false, 0);


	/******************** Bremsttarget ******************/
	G4Box *Bremstarget_solid = new G4Box("Bremstarget_solid", bremstarget_edge_length * 0.5, bremstarget_edge_length * 0.5, bremstarget_thickness * 0.5);

	G4LogicalVolume *Bremstarget_logical = new G4LogicalVolume(Bremstarget_solid, gold, "Bremstarget_logical", 0, 0, 0);

	//Visualisierung (Farbe)
	Bremstarget_logical->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));

	G4VPhysicalVolume *Bremstarget_physical = new G4PVPlacement(0, G4ThreeVector(0, 0, -detector_to_bremstarget/2), Bremstarget_logical, "Bremstarget", World_logical, false, 0);


	/******************** Detector ******************/
	G4Tubs *Detector_solid = new G4Tubs("Detector_solid", 0, detector_radius, detector_length * 0.5, 0, twopi);

	G4LogicalVolume *Detector_logical = new G4LogicalVolume(Detector_solid, detector_material, "Detector_logical", 0, 0, 0);

	//Visualisierung (Farbe)
	Detector_logical->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));

	G4VPhysicalVolume *Detector_physical = new G4PVPlacement(0, G4ThreeVector(0, 0, detector_to_bremstarget/2), Detector_logical, "Detector", World_logical, false, 0);

	print_info();
	return World_physical;
}

// Definiere das Detektorvolumen als Detektor/sensitives Volumen in Geant4
void DetectorConstruction::ConstructSDandField() {

	// TODO: Ask UFG whether this is the right detector, does this account for multiple gamma hits in the same event?
	// Do we get two entries in the ROOT file for the two gammas in the same event or are they summed / one is ignored?
	EnergyDepositionSD *DetectorSD = new EnergyDepositionSD("Detector_logical", "Detector_logical");
	G4SDManager::GetSDMpointer()->AddNewDetector(DetectorSD);
	DetectorSD->SetDetectorID(0);
	SetSensitiveDetector("Detector_logical", DetectorSD, true);
}

void DetectorConstruction::print_info() const
{
	printf("==============================================================\n");
	printf("  DetectorConstruction: Info (all dimensions in mm)\n");
	G4cout << G4endl;
	printf("> World dimensions             : ( %5.2f, %5.2f, %5.2f )\n", World_x, World_y, World_z);
	printf("==============================================================\n");
}
