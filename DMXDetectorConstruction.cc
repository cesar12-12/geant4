// --------------------------------------------------------------
//   GEANT 4 - Germanium Detector at IFUNAM
//									
//      For information related to this code contact: 
//	José Luis Hernández and Osmany González Reina
//      e-mail: luis.hernandez@ciencias.unam.mx
// --------------------------------------------------------------
// Comments
//
//                     Germanium Detector
//                    (17th September 2018)
//
// DetectorConstruction program
// --------------------------------------------------------------

#include "DMXDetectorConstruction.hh"
#include "DMXDetectorMessenger.hh"
#include "DMXScintSD.hh"
#include "DMXPmtSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Paraboloid.hh"

#include <math.h>
#define PI 3.14159265

DMXDetectorConstruction::DMXDetectorConstruction()  
{


  theUserLimitsForRoom     = 0; 
  theUserLimitsForDetector = 0; 
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theRoomTimeCut      = 1000. * nanosecond;
  theMinEkine         = 250.0*eV;
  theRoomMinEkine     = 250.0*eV;
  LXeSD.Put(0);                   

}

DMXDetectorConstruction::~DMXDetectorConstruction() 
{
  delete theUserLimitsForRoom;
  delete theUserLimitsForDetector;
  //delete detectorMessenger;
}

void DMXDetectorConstruction::DefineMaterials() 
{

 #include "DMXDetectorMaterial.icc"

}

G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  DefineMaterials();
 
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.5, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.9, 0.2, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;


  G4double wallThick   = 24.*cm;
  G4double worldWidth  = 600.0*cm + 2.*wallThick; 
  G4double worldLength = 600.0*cm + 2.*wallThick; 
  G4double worldHeight = 600.0*cm + 2.*wallThick; 
 
  G4bool checkOverlaps = true;

  G4RotationMatrix* rotation = new G4RotationMatrix();
  rotation->rotateX(90.*deg);
  rotation->rotateY(0.*deg);
  rotation->rotateZ(0.*deg);

  G4Box* world_box = new G4Box("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  G4double worldR = 115*cm;
  G4double labR =  worldR-10*cm;

  // G4Sphere* world_box = new G4Sphere("world_box",0.0*cm,worldR,0.0*deg,360*deg,0.0*deg,180*deg);
  world_log  = new G4LogicalVolume(world_box, world_mat, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"world_phys", world_log, NULL, false, checkOverlaps);
  G4VisAttributes* world_vat= new G4VisAttributes(white);
  world_vat->SetVisibility(true);
  world_log->SetVisAttributes(world_vat);

  G4double labWidth  = worldWidth  - 2.*wallThick; 
  G4double labLength = worldLength - 2.*wallThick;
  G4double labHeight = worldHeight - 2.*wallThick; 

  G4Box* lab_box = new G4Box("lab_box", 0.5*labWidth, 0.5*labLength, 0.5*labHeight );
  //G4Sphere* lab_box =new G4Sphere("lab_box",0.0*cm,labR,0.0*deg,360*deg,0.0*deg,180*deg);
  lab_log  = new G4LogicalVolume(lab_box, lab_mat, "lab_log");
  lab_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "lab_phys",lab_log, world_phys, false,checkOverlaps);
  G4VisAttributes* lab_vat= new G4VisAttributes(white);
  lab_vat->SetVisibility(true);
  lab_log->SetVisAttributes(lab_vat);

  //----------------------
  //------------------------------                                                                                     
  // HPGe Detector                                                                                                     
  //------------------------------                                                                                     
  
  
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  // General and scale parameters                                                                                      
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  
  G4double phimin = 0.0*deg; // minumum  angle                                                                   
  G4double phimax = 360*deg; // Maximum angle (complete cylinder)                                                
  G4double beta = 2.55/2.52; // scale factor (coin in radiography)                                               
  G4double alpha = 5.13/5.1; // scale factor (drawn geometry)                                                    
  G4double rmin1 = 0.0*cm; // Inner standar raius                                                                
  G4double rvacc = 0.1*cm; // gap between PopTop and lead shield
  
  G4double rminHPGe = 0.78*cm; // Borus contact inner radius Nominal rminHPGe = 0.48*cm;                                                     
  G4double rmaxHPGe = alpha*0.5*4.99*cm; // HPGe maximum radius Nominal alpha*0.5*5.09*cm                                                  
  G4double heightHPGe = 0.25*alpha*5.10*cm; // HPGe height Nominal 0.25*alpha*5.13*cm                                                       
  G4double HeightHPGe = 4*heightHPGe; // HPGe height                                                             
  G4double DLwidth = 0.65*mm; // Width of Dead layer    Nominal = 0.65*mm Max = 0.75*mm
  G4double DLwidth2 = 1.50*mm; //dead layer for the "side" part of the crystal
  G4double rhoDifHPGe = 0.7*cm;
  G4double heightHoleHPGe = 0.85*HeightHPGe;
  G4double sigHPGe = 0.5*HeightHPGe-rhoDifHPGe-DLwidth;
  G4double heightHolePos = 0.5*heightHoleHPGe - 2*heightHPGe;
  G4double rDLout = rmaxHPGe;
  G4double rDLin = rmaxHPGe - DLwidth;
  G4double DLheight = 0.5*(HeightHPGe - rhoDifHPGe - DLwidth);
  G4double DLpos1 = sigHPGe - DLheight;
  G4int nHPGe_3 = 50;
  G4double deltaDL = (rhoDifHPGe+DLwidth)/nHPGe_3;
  double zHPGe_3[nHPGe_3+1], rOutHPGe_3[nHPGe_3+1], rInHPGe_3[nHPGe_3+1];

  for(int i=0 ; i<=nHPGe_3; i++){
    
    zHPGe_3[i] = sigHPGe + i*deltaDL;
    rOutHPGe_3[i] = rmaxHPGe-rhoDifHPGe-DLwidth+sqrt(pow(rhoDifHPGe+DLwidth,2)-pow(zHPGe_3[i]-sigHPGe,2));

    if(zHPGe_3[i]<=0.5*HeightHPGe-DLwidth)
      //rInHPGe_3[i] = rmaxHPGe-rhoDifHPGe-DLwidth+sqrt(pow(rhoDifHPGe,2)-pow(zHPGe_3[i]-sigHPGe,2));
      rInHPGe_3[i] = 0.*cm;
    else
      rInHPGe_3[i] = 0.*cm;
  };

  G4int nHPGe_2 = 50;
  G4double deltaHPGe = rhoDifHPGe/nHPGe_2;
  G4double zHPGe_2[nHPGe_2+1], rInHPGe_2[nHPGe_2+1], rOutHPGe_2[nHPGe_2+1];
  G4int nHPGe_1 = 4;
  G4double zHPGe_1[nHPGe_1], rInHPGe_1[nHPGe_1], rOutHPGe_1[nHPGe_1];


  if(heightHoleHPGe<HeightHPGe-rhoDifHPGe-DLwidth){
    zHPGe_1[0] = -0.5*HeightHPGe;
    zHPGe_1[1] = -0.5*HeightHPGe+heightHoleHPGe;
    zHPGe_1[2] = -0.5*HeightHPGe+heightHoleHPGe;
    zHPGe_1[3] = sigHPGe;
    rInHPGe_1[0] = rminHPGe;
    rInHPGe_1[1] = rminHPGe;
    rInHPGe_1[2] = 0*cm;
    rInHPGe_1[3] = 0*cm;
    for(int j=0; j<nHPGe_1 ; j++)
      rOutHPGe_1[j] = rmaxHPGe-DLwidth2;
    for(int i=0; i<=nHPGe_2; i++){
      zHPGe_2[i] = sigHPGe +i*deltaHPGe;
      rOutHPGe_2[i] = rmaxHPGe-rhoDifHPGe-DLwidth+sqrt(pow(rhoDifHPGe,2)-pow(zHPGe_2[i]-sigHPGe,2));
      rInHPGe_2[i] = 0*cm;
    }
  }
  else{
    zHPGe_1[0] = -0.5*HeightHPGe;
    zHPGe_1[1] = -0.25*HeightHPGe;
    zHPGe_1[2] = 0*cm;
    zHPGe_1[3] = sigHPGe;
    for(int j=0; j<nHPGe_1 ; j++){
      rInHPGe_1[j]= rminHPGe;
      rOutHPGe_1[j] = rmaxHPGe-DLwidth2;
    }

    for(int i=0; i<=nHPGe_2; i++){
      zHPGe_2[i] = sigHPGe +i*deltaHPGe;
      rOutHPGe_2[i] = rmaxHPGe-rhoDifHPGe-DLwidth+sqrt(pow(rhoDifHPGe,2)-pow(zHPGe_2[i]-sigHPGe,2));
      if(zHPGe_2[i]<=heightHoleHPGe-0.5*HeightHPGe)
	rInHPGe_2[i] = rminHPGe;
      else
	rInHPGe_2[i] = 0*cm;
    }
  }


  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  // Inner and outer contacts parameters                                                                               
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   

  G4double BorContactWidth = beta*0.01*cm; // Borus inner contact width                                          
  G4double LiContactWidth = beta*0.02*cm; // Ltihium outer contact width                                         
  G4double LiContactHeight = 0.5*10.92*cm;
  G4double LiContactPos = 0.5*HeightHPGe-LiContactHeight;
  G4double LiContactPos2 = 0.5*HeightHPGe+0.5*LiContactWidth;
  G4double LiContactPos3 = LiContactPos2 - 2*LiContactHeight - LiContactWidth;
  G4double EndGapWidth = beta*0.1*cm; // EndGap Width                                                            
  G4double posEndGap = -2*heightHPGe-EndGapWidth; // EndGap inferior position                                    
  G4double AlWindowra = 4*rminHPGe; // Radius of aluminium window                                                
  G4double AlWindowhe = 0.5*EndGapWidth; // Window height                                                        

  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  // Cool finger parameters                                                                                            
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   

  G4double Coldheight = beta*2.5*cm; // Height of cool finger                                                    

  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  // Cryostat parameters                                                                                               
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   

  G4double rAlCover = 0.5*6.8*cm; // Radius of Aluminium cover                                                   
  G4double AlCoverWidth = beta*0.5*mm; // Width of Aluminium cover                                               
  G4double AlCoverHeight = beta*0.5*14.23*cm; // Height of Aluminium cover                                       
  G4double AlCoverPos = (rAlCover-rmaxHPGe)+2*heightHPGe - AlCoverHeight; // Aluminum cover position             
  G4double AlTopPos = AlCoverPos+AlCoverHeight;
  G4double deltaAux = AlTopPos - AlCoverWidth - LiContactPos2- 0.5*LiContactWidth;

  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  // Lead shield parameters                                                                                            
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   

  G4double rminShield1 = 7.75*cm; // Inner radius of lead shield                                                 
  G4double Shield1Width = 4.75*cm; // Width of lead shield                                                       
  G4double rmaxShield1 = rminShield1 + Shield1Width; // Maximum radius of lead shield                            
  G4double Topheight = 21*cm;
  G4double heightShield1 = 4*12.7*cm; // Height of Pb shield                                                     
  G4double Shield1pos = AlTopPos - Topheight + 0.5*heightShield1; // Pb shield position                          
  G4double rminShield2 = rmaxShield1; // Inner radius of Cu shield                                               
  G4double Shield2Width1 = 1.9*mm; // Width of Cu shield                                                         
  G4double Shield2Width2 = 1.4*mm; // Mid-shield width                                                           
  G4double rmidShield2 = rminShield2 + Shield2Width1;//                                                          
  G4double rmaxShield2 = rmidShield2 + Shield2Width2; //                                                         
  G4double rminShield3 = rmaxShield2 + rvacc; //                                                                 
  G4double Shield3Width = 6.6*mm; //                                                                             
  G4double rmaxShield3 = rminShield3 + Shield3Width; //                                                          
  G4double HeadShieldWidth = 6*cm;
  G4double HeadPos = Shield1pos + 0.5*heightShield1 + 0.5*HeadShieldWidth; //  
  
  G4double rminCu = 38*mm; //Min radius lining Cu part 1
  G4double LinCuwidth = 1.6*mm; //Width of Cu Graded Lining
  G4double rlinCu= rminCu + LinCuwidth; //Max radius of Graded Lining Cu
  G4double ShieldPbwidth1 = 204.9*mm; //Width of Pb Part 1
  G4double rPbShield = rlinCu + ShieldPbwidth1; // Max radius Shield Pb Part 1
  G4double ShieldStwidth = 9.5*mm; // Width of Outer Jacket Carbon Steel
  G4double rStShield = rPbShield + ShieldStwidth; // Max Radius Steel Jacket Part 1
  G4double rminCu2 = 139.5*mm; //Min radius lining Cu part2
  G4double rlinCu2 = rminCu2 + LinCuwidth; //Max radius of Graded Lining Cu 2
  G4double ShieldPbwidth2 = 103.4 mm; //Width of Pb Part 2
  G4double rPbshield2 = rlinCu2 + ShieldPbwidth2; //Max radius Shield Pd Part 2
  G4double rStShield2 = rPbshield2 + ShieldStwidth; // Max radius Steel Jacket Part 2
  G4double heightpart1 = 0.5*116*mm; // Height Part 1 Shield
  G4double heightpart2 = 0.5*407*mm; // Height Part 2 Shield
  G4double rmaxlead = 244.5*mm; // Max radius of Cu and Pb lead
  G4double heightleadPb = 101.9*mm; // Height of Pb lead
  G4double rmaxleadSt = rmaxlead + ShieldStwidth; // Max radius os Outer Jacket  Carbon Steel Lead


 // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  // Experimental setup parameters                                                                                     
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   

  G4double rSource = 0.75*cm; //                                                                                 
  G4double heightSource = 0.25*cm; //                                                                            
  G4double tableX = 30*cm;
  G4double tableY = 5.0*cm;
  G4double tablePos = AlTopPos - Topheight - 0.5*tableY;
  G4double legX = 0.5*tableY;
  G4double legY = 70.0*cm;
  G4double legPos = tablePos - 0.5*tableY - 0.5*legY;
  G4double rinDew1 = rAlCover+AlCoverWidth;
  G4double widthDew = 5.0*cm;
  G4double routDew1 = rinDew1+widthDew;
  G4double DewProp = 0.7;
  G4double Dewheight1 = (1-DewProp)*0.5*legY;
  G4double DewPos1 = tablePos - 0.5*tableY - Dewheight1;
  G4double rinDew2 = 15.0*cm;
  G4double routDew2 = rinDew2 + widthDew;
  G4double Dewheight2 = DewProp*0.5*legY-widthDew;
  G4double DewPos2 = DewPos1-Dewheight1-Dewheight2;
  G4double DewPos3 = DewPos2-Dewheight2-0.5*widthDew;

  // Electronics complements                                                                                           
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   


  G4double aux1Elec = AlTopPos - 2*AlCoverWidth;
  G4double aux2Elec = -(tablePos+0.5*tableY);
  G4double heightElec = 0.5*(aux2Elec - (2*AlCoverHeight - aux1Elec));
  G4double radElec = rAlCover + 0.5*cm;
  G4double discElecWidth = 0.5*0.5*cm;
  G4double discElecRad = 1.2*cm;
  G4double heightElecFin = 0.5*(HeightHPGe-2*EndGapWidth-2*discElecWidth);
  G4double radElecFin = 0.5*0.77*cm;
  G4double zPosElecFin = posEndGap - EndGapWidth - heightElecFin;
  G4double zPosElecDisc = zPosElecFin - heightElecFin - discElecWidth;
  G4double zPosElec = AlCoverPos - AlCoverHeight - heightElec;
  G4double heightCool = 2*AlCoverHeight - AlCoverWidth - deltaAux - LiContactWidth - HeightHPGe - 2*EndGapWidth - 2*heightElecFin - 2*discElecWidth + 2*heightElec;
  G4double posCool = zPosElecDisc - discElecWidth - 0.5*heightCool;
  G4double rCool = 0.5*1.5*cm; // Radius of cool finger           

  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                                          
  // Sample Container                                                                                                                         
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                                          

  G4double SCheight = 4.2*cm;
  G4double SCwidth = 0.9*mm;
  G4double SCrad_out = 0.5*4.75*cm;
  G4double SCrad_in = SCrad_out - SCwidth;
  G4double SCdist = 0.0*cm;
  G4double SCpos1 = AlTopPos + SCdist + 0.5*SCwidth;
  G4double SCpos2 = AlTopPos + 0.5*SCheight;
  G4double SCpos3 = SCpos2 + 0.5*SCheight-0.5*SCwidth;
  bool opSC = false;
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                                          
  // Sample Source                                                                                                                            
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                                          

  G4double Sampleheight = SCheight - 2*SCwidth;


  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   
  // SENSIBLE DETECTOR CONSTRUCTION                                                                                    
  // -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-                                                                                   

  //Aluminum Cover
  //
  // G4Tubs* coolfingerhole = new G4Tubs("coolfingerhole",0,rCool+0.01*mm,0.5*heightCool,phimin,phimax);

  G4Tubs* solidAlCover = new G4Tubs("solidAlCover",0,rAlCover+AlCoverWidth,AlCoverHeight,phimin,phimax);
  //G4SubtractionSolid *solidAlCover = new G4SubtractionSolid("solidAlCover",
  //						    presolidAlCover,coolfingerhole,0,G4ThreeVector(0.,0.,-0.5*AlCoverHeight));
  logicAlCover = new G4LogicalVolume(solidAlCover, panel_mat, "logicAlCover");
  physiAlCover = new G4PVPlacement(0, G4ThreeVector(0.,0.,AlCoverPos), "physiAlCover", logicAlCover, lab_phys, false, checkOverlaps);
  

  //Vacuum inside cover
  G4Tubs* solidAlCoverVac = new G4Tubs("solidAlCoverVac",0,rAlCover,AlCoverHeight-AlCoverWidth,phimin,phimax);
  logicAlCoverVac = new G4LogicalVolume(solidAlCoverVac, vacuum_mat, "logicAlCoverVac");
  physiAlCoverVac = new G4PVPlacement(0, G4ThreeVector(0.,0.,0), "physiAlCoverVac", logicAlCoverVac, physiAlCover, false, checkOverlaps);
  //Lithium contact
  G4double lithContactPos = 0.88467290*cm;//-AlCoverPos-LiContactPos;//-3.7795847*cm;
  G4Tubs * solidLiContact = new G4Tubs("solidLiContact",0,rmaxHPGe+LiContactWidth,LiContactHeight+LiContactWidth,phimin,phimax);
  logicLiContact = new G4LogicalVolume( solidLiContact, LiContactMater, "logicLiContact");
  physiLiContact = new G4PVPlacement(0, G4ThreeVector(0.,0.,lithContactPos), 
				     "physiLiContact", logicLiContact, physiAlCoverVac, false, checkOverlaps);
  G4VisAttributes* LiContactcolor = new G4VisAttributes(yellow); 
  LiContactcolor->SetVisibility(true);
  logicLiContact->SetVisAttributes(LiContactcolor);
  //Vac inside contact
  G4Tubs * solidLiContactVac = new G4Tubs("solidLiContactVac",0,rmaxHPGe,LiContactHeight,phimin,phimax);
  logicLiContactVac = new G4LogicalVolume( solidLiContactVac, vacuum_mat, "logicLiContactVac");
  physiLiContactVac = new G4PVPlacement(0, G4ThreeVector(0.,0.,0), 
				     "physiLiContactVac", logicLiContactVac, physiLiContact, false, checkOverlaps);
  //Dead layer cylindrical section of Ge crystal
  G4double dlpos1 =-LiContactPos+DLpos1;
  G4Tubs* solidDlayer1 = new G4Tubs("solidDlayer1",0, rDLout, DLheight, phimin, phimax);
  logicDlayer1 = new G4LogicalVolume( solidDlayer1, BEGe_mat, "logicDlayer1");
  physiDlayer1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,dlpos1), 
				   "physiDlayer1", logicDlayer1, physiLiContactVac, false, checkOverlaps);
  //cylindrical section of Ge crystal
  G4Polycone* solidDisc1 = new G4Polycone("solidDisc1",phimin, phimax, nHPGe_1, zHPGe_1, rInHPGe_1, rOutHPGe_1);
  logicDisc1 = new G4LogicalVolume(solidDisc1, BEGe_mat, "logicDisc1");
  physiDisc1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,-DLpos1), 
				 "physiDisc1", logicDisc1, physiDlayer1, false, checkOverlaps);
  //Dead layer of round section of Ge crystal
  G4Polycone* solidDlayer2 = new G4Polycone("solidDlayer2", phimin, phimax, nHPGe_3, zHPGe_3, rInHPGe_3, rOutHPGe_3);
  logicDlayer2 = new G4LogicalVolume( solidDlayer2, BEGe_mat, "logicDlayer2");
  physiDlayer2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,-LiContactPos), 
				   "physiDlayer2", logicDlayer2, physiLiContactVac, false, checkOverlaps);
  //Round section of Ge crystal
  G4Polycone* solidDisc2 = new G4Polycone("solidDisc2",phimin, phimax, nHPGe_2, zHPGe_2, rInHPGe_2, rOutHPGe_2);
  logicDisc2 = new G4LogicalVolume( solidDisc2, BEGe_mat, "logicDisc2");
  physiDisc2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "physiDisc2",
				 logicDisc2, physiDlayer2, false, checkOverlaps);
  G4VisAttributes* solidDisccolor = new G4VisAttributes(red);
  solidDisccolor->SetVisibility(true);
  logicDisc2->SetVisAttributes(solidDisccolor);
  logicDisc1->SetVisAttributes(solidDisccolor);

  //Boron contact
  G4double contactPos = -DLpos1+heightHolePos;
  G4Tubs* solidContact = new G4Tubs("solidContact", 0, rminHPGe, 0.5*heightHoleHPGe, phimin, phimax);//rmin=rminHPGe-BorContactWidth
  logicContact = new G4LogicalVolume(solidContact, ContactMater, "logicContact");
  physiContact = new G4PVPlacement(0, G4ThreeVector(0.,0.,contactPos), 
				   "physiContact", logicContact, physiDlayer1, false, checkOverlaps);
  G4VisAttributes* contactcolor = new G4VisAttributes(magenta);
  contactcolor->SetVisibility(true);
  logicContact->SetVisAttributes(contactcolor);
  //Vac inside Boron Contact
  G4Tubs* solidContactVac = new G4Tubs("solidContactVac",0,rminHPGe-BorContactWidth, 0.5*heightHoleHPGe, phimin, phimax);
  logicContactVac = new G4LogicalVolume(solidContactVac,vacuum_mat , "logicContactVac");
  physiContactVac = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
                                   "physiContactVac", logicContactVac, physiContact, false, checkOverlaps);

  //carbon endgap
  G4double endgappos = -LiContactPos + posEndGap;
  G4Tubs* solidEndgap = new G4Tubs("solidEndgap",rmin1, rmaxHPGe, EndGapWidth, phimin, phimax);
  logicEndgap = new G4LogicalVolume(solidEndgap, EndgapMater, "logicEndgap");
  physiEndgap = new G4PVPlacement(0, G4ThreeVector(0.,0.,endgappos), 
				  "physiEndgap", logicEndgap, physiLiContactVac, false, checkOverlaps);
  //Cooling finger
  G4double elecfingerpos = -LiContactPos+zPosElecFin;
  G4Tubs* solidEFinger = new G4Tubs("solidEFinger",rmin1,radElecFin,heightElecFin,phimin,phimax);
  logicEFinger = new G4LogicalVolume( solidEFinger, panel_mat, "logicEFinger");
  physiEFinger = new G4PVPlacement(0, G4ThreeVector(0.,0.,elecfingerpos), "physiEFinger", logicEFinger,physiLiContactVac, false, checkOverlaps);

  G4double elecdiscpos = -LiContactPos+zPosElecDisc;
  G4Tubs* solidEDisc = new G4Tubs("solidEDisc",rmin1,discElecRad,discElecWidth,phimin,phimax);
  logicEDisc = new G4LogicalVolume(solidEDisc, panel_mat, "logicEDisc", 0,0,0);
  physiEDisc = new G4PVPlacement(0, G4ThreeVector(0.,0.,elecdiscpos),"physiEDisc",logicEDisc,physiLiContactVac, false, checkOverlaps);
  
  G4double coolfinheigh1 = LiContactHeight + elecdiscpos - discElecWidth;
  G4Tubs* coolfin1 = new G4Tubs("coolfin1",0,rCool,0.5*coolfinheigh1,phimin,phimax);
  logicCoolFin1 = new G4LogicalVolume(coolfin1, panel_mat, "logicCoolFin1");
  G4double coolfin1Pos = - LiContactHeight + 0.5*coolfinheigh1; 
  physiCoolFin1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,coolfin1Pos),"physiCoolFin1",logicCoolFin1,physiLiContactVac, false, checkOverlaps);

  G4double coolfinheigh2 = LiContactWidth;
  G4Tubs* coolfin2 = new G4Tubs("coolfin2",0,rCool,0.5*coolfinheigh2,phimin,phimax);
  logicCoolFin2 = new G4LogicalVolume(coolfin2, panel_mat, "logicCoolFin2");
  G4double coolfin2Pos = - LiContactHeight-LiContactWidth + 0.5*coolfinheigh2;
  physiCoolFin2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,coolfin2Pos),"physiCoolFin2",logicCoolFin2,physiLiContact, false, checkOverlaps);

  G4double coolfinheigh3 = AlCoverHeight-AlCoverWidth + lithContactPos - LiContactHeight-LiContactWidth;//+LiContactWidth;
  G4Tubs* coolfin3 = new G4Tubs("coolfin3",0,rCool,0.5*coolfinheigh3,phimin,phimax);
  logicCoolFin3 = new G4LogicalVolume(coolfin3, panel_mat, "logicCoolFin3");
  G4double coolfin3Pos = -AlCoverHeight+AlCoverWidth + 0.5*coolfinheigh3;// -0.5*LiContactWidth;
  physiCoolFin3 = new G4PVPlacement(0, G4ThreeVector(0.,0.,coolfin3Pos),"physiCoolFin3",logicCoolFin3,physiAlCoverVac, false, checkOverlaps);

  G4double coolfinheigh4 = AlCoverWidth;
  G4Tubs* coolfin4 = new G4Tubs("coolfin4",0,rCool,0.5*coolfinheigh4,phimin,phimax);
  logicCoolFin4 = new G4LogicalVolume(coolfin4, panel_mat, "logicCoolFin4");
  G4double coolfin4Pos = -AlCoverHeight+0.5*coolfinheigh4;
  physiCoolFin4 = new G4PVPlacement(0, G4ThreeVector(0.,0.,coolfin4Pos),"physiCoolFin4",logicCoolFin4,physiAlCover, false, checkOverlaps);

 
  G4Tubs* solidElectronics = new G4Tubs("solidElectronics",rCool,radElec,heightElec,phimin,phimax);
  logicElectronics = new G4LogicalVolume(solidElectronics, panel_mat, "logicElectronics");
  physiElectronics = new G4PVPlacement(0, G4ThreeVector(0.,0.,zPosElec),"physiElectronics",logicElectronics,lab_phys, false, checkOverlaps);
  
  G4Tubs* solidCool = new G4Tubs("cool",rmin1,rCool,heightElec,phimin,phimax);
  logicCool = new G4LogicalVolume(solidCool, panel_mat, "logicCool");
  physiCool = new G4PVPlacement(0,G4ThreeVector(0.,0.,zPosElec),"physiCool",logicCool,lab_phys, false, checkOverlaps);


  

  /* 
  //SAMPLE
  G4Tubs* solidSC1 = new G4Tubs("solidSC1",SCrad_in,SCrad_out,0.5*SCheight,phimin,phimax);
  logicSC1 = new G4LogicalVolume(solidSC1, pp_base_mat, "logicSC1");
  physiSC1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,SCpos2),"physiSC1",logicSC1, lab_phys,false, checkOverlaps);

  G4Tubs* solidSC2 = new G4Tubs("solidSC2",rmin1,SCrad_in,0.5*SCwidth,phimin,phimax);
  logicSC2 = new G4LogicalVolume( solidSC2, pp_base_mat, "logicSC2");
  physiSC2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,SCpos1),"physiSC2",logicSC2, lab_phys,false, checkOverlaps);

  physiSC3 = new G4PVPlacement(0, G4ThreeVector(0.,0.,SCpos3),"physiSC3",logicSC2, lab_phys,false, checkOverlaps);

  G4Tubs* solidSource = new G4Tubs("solidSource",rmin1,SCrad_in,0.5*Sampleheight,phimin,phimax);
  logicSource = new G4LogicalVolume( solidSource, pp_base_mat, "logicSource");
  physiSource = new G4PVPlacement(0,G4ThreeVector(0.,0.,SCpos2),"physiSource",logicSource,lab_phys,false, checkOverlaps);

  

  //Point source
  G4double sourceR = 12.75*mm;                                                                                                                    
  G4double sourceH = 3.2*mm;                
  G4RotationMatrix* sourceRot = new G4RotationMatrix;                                                                                             
  sourceRot->rotateY(PI/2.*rad);                                                                                                                
  //sourceRot->rotateZ(PI/2.*rad);
  //G4double sourcePosZ = pbBlockhole_placement+0.5*pbBlockholez+0.5*sourceH+0.01*mm;//                                                           
  //G4double sourcePosZ = 0.5*covercylinderheight + covertube_out_height + 0.01*mm + 8.95*cm - 0.01*mm - 0.5*sourceH + maxDist;                   
  G4Tubs* source = new G4Tubs("source", 0, sourceR, sourceH/2., 0.*deg, 360.*deg);                                                                
  //G4Box* source = new G4Box("source", sourcex/2 , sourcey/2 , sourcez/2);                                                                     
  logicSource  = new G4LogicalVolume(source, glass_port_mat, "logicSource");                                                                        
  physiSource = new G4PVPlacement(0, G4ThreeVector(0, 0, 28.42*cm),  "physiSource",logicSource, lab_phys, false,checkOverlaps);
  
  */
  //SHIELD
  G4Tubs* solidShield1 = new G4Tubs("solidShield1",rminShield1, rmaxShield1,0.5*heightShield1,phimin,phimax);
  logicShield1 = new G4LogicalVolume( solidShield1, Pbshield_mat, "logicShield1");
  physiShield1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,Shield1pos),"physiShield1",logicShield1,lab_phys,false, checkOverlaps);

  G4Tubs* solidShield2 = new G4Tubs("solidShield2",rminShield2, rmidShield2,0.5*heightShield1,phimin,phimax);
  logicShield2 = new G4LogicalVolume( solidShield2, Cushield_mat, "logicShield2");
  physiShield2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,Shield1pos),"physiShield2",logicShield2,lab_phys,false, checkOverlaps);

  G4Tubs* solidShield21 = new G4Tubs("solidShield21",rmidShield2, rmaxShield2,0.5*heightShield1,phimin,phimax);
  logicShield21 = new G4LogicalVolume(solidShield21, shieldcover_mat, "logicShield21");
  physiShield21 = new G4PVPlacement(0, G4ThreeVector(0.,0.,Shield1pos),"physiShield21",logicShield21,lab_phys,false, checkOverlaps);

  G4Tubs* solidShield3 = new G4Tubs("solidShield21",rminShield3, rmaxShield3,0.5*heightShield1,phimin,phimax);
  logicShield3 = new G4LogicalVolume(solidShield3, shieldcover_mat, "logicShield3");
  physiShield3 = new G4PVPlacement(0, G4ThreeVector(0.,0.,Shield1pos),"physiShield3", logicShield3,lab_phys,false, checkOverlaps);
  
  G4Tubs* solidShieldH = new G4Tubs("solidShieldH",rmin1, rmaxShield3,0.5*HeadShieldWidth,phimin,phimax);
  logicShieldH = new G4LogicalVolume(solidShieldH, shieldcover_mat, "logicShieldH");
  physiShieldH = new G4PVPlacement(0, G4ThreeVector(0.,0.,HeadPos),"physiShieldH",logicShieldH,lab_phys,false, checkOverlaps);

  //NEW SHIELD CONSTRUCTION

  G4Tubs* LiningCu1 = new G4Tubs("LiningCu1", rminCu, rlinCu, heightpart1, phimin, phimax);
  logicLin1 = new G$LogicalVolume( LinCu1, Cushield_mat, "logicLin1");  
  //TABLE

  G4Box* solidTable = new G4Box("solidTable",tableX,tableX,0.5*tableY);
  logicTable = new G4LogicalVolume(solidTable, panel_mat, "logicTable");
  physiTable = new G4PVPlacement(0, G4ThreeVector(0.,0.,tablePos),"physiTable",logicTable,lab_phys,false, checkOverlaps);
  G4Box* tableHole = new G4Box("tableHole", 6.5*cm, 6.5*cm, 0.5*tableY );
  tableHole_log = new G4LogicalVolume(tableHole, lab_mat, "tableHole_log");
  tableHole_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"tableHole_phys",tableHole_log,physiTable,false, checkOverlaps);
  G4Tubs* tableFinger = new G4Tubs("tableFinger",0.*cm,1.59*cm,0.5*tableY, phimin,phimax);
  tableFinger_log = new G4LogicalVolume(tableFinger, panel_mat, "tableFinger_log");
  tableFinger_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"tableFinger_phys",tableFinger_log,tableHole_phys,false, checkOverlaps);
  
  G4Box* solidLeg = new G4Box("solidLeg",legX,legX,0.5*legY);
  logicLeg = new G4LogicalVolume(solidLeg, panel_mat, "logicLeg");
  physiLeg1 = new G4PVPlacement(0, G4ThreeVector(tableX-legX,tableX-legX,legPos),
				"physiLeg1",logicLeg,lab_phys,false, checkOverlaps);
  physiLeg2 = new G4PVPlacement(0, G4ThreeVector(-tableX+legX,tableX-legX,legPos),
				"physiLeg2",logicLeg,lab_phys,false, checkOverlaps);
  physiLeg3 = new G4PVPlacement(0, G4ThreeVector(-tableX+legX,-tableX+legX,legPos),
				"physiLeg3",logicLeg,lab_phys,false, checkOverlaps);
  physiLeg4 = new G4PVPlacement(0, G4ThreeVector(tableX-legX,-tableX+legX,legPos),
				"physiLeg4",logicLeg,lab_phys,false, checkOverlaps);

  /*
  G4Cons* solidDewar1 = new G4Cons("solidDewar1",rinDew2,routDew2,rinDew1,routDew1,Dewheight1,phimin,phimax);
  logicDewar1 = new G4LogicalVolume( solidDewar1, panel_mat, "logicDewar1");
  physiDewar1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,DewPos1),
				  "physiDewar1",logicDewar1,lab_phys,false, checkOverlaps);

  G4Tubs* solidDewar2 = new G4Tubs("solidDewar2",rinDew2,routDew2,Dewheight2,phimin,phimax);
  logicDewar2 = new G4LogicalVolume( solidDewar2, panel_mat, "logicDewar2");
  physiDewar2 = new G4PVPlacement(0,G4ThreeVector(0.,0.,DewPos2),
				  "physiDewar2",logicDewar2,lab_phys,false, checkOverlaps);

  G4Tubs* solidDewar3 = new G4Tubs("solidDewar3",rmin1,routDew2,0.5*widthDew,phimin,phimax);
  logicDewar3 = new G4LogicalVolume( solidDewar3, panel_mat, "logicDewar3");
  physiDewar3 = new G4PVPlacement(0,G4ThreeVector(0.,0.,DewPos3),
				  "physiDewar3",logicDewar3,lab_phys,false, checkOverlaps);
  */


  const G4int dewar_num = 3;
  G4double dewar_rmax[dewar_num] = {routDew2, routDew2, routDew1};
  G4double dewar_rmin[dewar_num] = {0.*cm,0.*cm,0.*cm};
  G4double dewarOutheight1 = widthDew+Dewheight2*2;
  G4double dewarOutTotalheight = 2*Dewheight1+widthDew+Dewheight2*2;
  G4double dewar_z[dewar_num] = {0.*cm,dewarOutheight1,dewarOutTotalheight};
  G4Polycone* dewar = new G4Polycone("dewar", 0.*deg, 360.*deg, dewar_num, dewar_z, dewar_rmin, dewar_rmax);
  dewar_log = new G4LogicalVolume( dewar, panel_mat, "dewar");
  G4double dewarPos = DewPos3-(0.5*widthDew);
  dewar_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,dewarPos),
				 "dewar_phys",dewar_log,lab_phys,false, checkOverlaps);

  const G4int dewarair_num = 5;
  G4double dewarair_rmax[dewarair_num] = {routDew2-3.*mm, routDew2-3.*mm, routDew1-3.*mm, rinDew1,rinDew1};
  G4double dewarair_rmin[dewarair_num] = {0.*cm,0.*cm,0.*cm,0.*cm,0.*cm,};
  G4double dewarair_z[dewarair_num] = {0.*cm,dewarOutheight1-3.*mm,dewarOutTotalheight-6.*mm,dewarOutTotalheight-6.*mm,dewarOutTotalheight-3.*mm};
  G4Polycone* dewarair = new G4Polycone("dewarair", 0.*deg, 360.*deg, dewarair_num, dewarair_z, dewarair_rmin, dewarair_rmax);
  dewarair_log = new G4LogicalVolume(dewarair,vacuum_mat, "dewarair");
  G4double dewarairPos = 3.*mm;
  dewarair_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,dewarairPos),
                                 "dewarair_phys",dewarair_log,dewar_phys,false, checkOverlaps);

  
  const G4int dewarin_num = 4;
  G4double dewarin_rmax[dewarin_num] = {rinDew2, rinDew2, rinDew1, rinDew1};
  G4double dewarin_rmin[dewarin_num] = {0.*cm,0.*cm,0.*cm, 0.*cm};
  G4double dewarin_z[dewarin_num] = {0.*cm,dewarOutheight1-2*widthDew-6.*mm,dewarOutheight1-2*widthDew-6.*mm,dewarOutTotalheight-widthDew-3.*mm};
  G4Polycone* dewarin = new G4Polycone("dewarin", 0.*deg, 360.*deg, dewarin_num, dewarin_z, dewarin_rmin, dewarin_rmax);
  dewarin_log = new G4LogicalVolume(dewarin, panel_mat, "dewarin_log");
  G4double dewarinPos = widthDew;
  dewarin_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,dewarinPos),
				    "dewarin_phys",dewarin_log,dewarair_phys,false, checkOverlaps);

  G4double dewarinLN_rmax[dewarin_num] = {rinDew2-3.*mm, rinDew2-3.*mm, rinDew1-3.*mm, rinDew1-3.*mm};
  G4double dewarinLN_rmin[dewarin_num] = {0.*cm,0.*cm,0.*cm, 0.*cm};
  G4double dewarinLN_z[dewarin_num] = {0.*cm,dewarOutheight1-2*widthDew-12.*mm,dewarOutheight1-2*widthDew-12.*mm,dewarOutTotalheight-widthDew-6.*mm};
  G4Polycone* dewarinLN = new G4Polycone("dewarinLN", 0.*deg, 360.*deg, dewarin_num, dewarinLN_z, dewarinLN_rmin, dewarinLN_rmax);
  dewarinLN_log = new G4LogicalVolume(dewarinLN, LN2_mat, "dewarinLN_log");
  G4double dewarinPosLN = 3.*mm;
  dewarinLN_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,dewarinPosLN),
				   "dewarinLN_phys",dewarinLN_log,dewarin_phys,false, checkOverlaps);
  G4VisAttributes* dewarLNcolor = new G4VisAttributes(blue);
  dewarLNcolor->SetVisibility(true);
  dewarinLN_log->SetVisAttributes(dewarLNcolor);

  G4Tubs* dewarFinger = new G4Tubs("dewarFinger", 0.*cm,1.59*cm,0.5*51.0*cm, phimin,phimax);
  dewarFinger_log = new G4LogicalVolume(dewarFinger, panel_mat, "dewarFinger_log");
  G4double dewarFinPos = dewarOutTotalheight-widthDew-6.*mm - 0.5*51.0*cm;
  dewarFinger_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,dewarFinPos),
				       "dewarFinger_phys",dewarFinger_log,dewarinLN_phys,false, checkOverlaps);

  return world_phys;

}

void DMXDetectorConstruction::ConstructSDandField()	
{
    if (LXeSD.Get() == 0)
      {
        G4String name="/DMXDet/LXeSD";
        DMXScintSD* aSD = new DMXScintSD(name);
        LXeSD.Put(aSD);
      }
    G4SDManager::GetSDMpointer()->AddNewDetector(LXeSD.Get());
    //    if (logicDisc1)
    SetSensitiveDetector(logicDisc1,LXeSD.Get());
      //if (logicDisc2)
    SetSensitiveDetector(logicDisc2,LXeSD.Get());
 
    return;

}




