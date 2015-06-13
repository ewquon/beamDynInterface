/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "beamDynInterfacePointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "displacementMotionSolver.H"

//#include "beamDyn.H" // BD namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (2)" << endl;
}


beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (3)" << endl;
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const beamDynInterfacePointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (4)" << endl;
}


beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const beamDynInterfacePointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (2,ptf)" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void beamDynInterfacePointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //const labelList& meshPoints = patch().meshPoints(); // returns polyPatch_.meshPoints(), i.e. node IDs
    const pointField& localPoints = patch().localPoints(); // returns polyPatch_.localPoints(), i.e. node coords

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();

// PARAMETERS FOR CAMBER/PLUNGE TEST
//    const vector amplitude(0,0.75,0);
//    const scalar omega(15.707963267948966);
//
//    // normalization guarantees that the max camber will be as specified
//    const scalar xoff(0.25);
//    const scalar chord(1);
//    //const vector camber(0,-0.25,0); // this blows up after ~0.06 seconds
//    //const vector camber(0,-0.2,0); // blows up after ~0.07 seconds
//    //const vector camber(0,-0.1,0); // this is convergent
//    const vector camber(0,-0.15,0); // this is convergent
//    scalar norm(xoff);
//    if( chord-xoff > xoff ) norm = chord-xoff;

    vectorList& disp = BD::disp();  // linear displacement
    vectorList& adisp = BD::adisp(); // angular displacement

    //Pout<< "Saved displacement : " << disp << endl;
    //All procs have the same value
    Info<< "Saved displacement : " << disp << endl;

    //double R[9];
    double ang;
    double tmp[3];
    vector v(vector::zero);
    vector a(vector::zero);

    //
    // --loop over all surface nodes
    //
    forAll(*this, ptI)
    {
//        /////////////////////////////////////////////////////////////////////////////////////////////
//        // TESTS:
//
//        // sinusoidal plunging motion, equivalent to Field<vector>::operator=(amplitude*sin(omega*t.value()));
//        //this->operator[](ptI) = amplitude*sin(omega*t.value());
//        
//        // sinusoidal plunging motion with variable camber
//        scalar xI(localPoints[ptI].component(0));
//        vector deform(camber);
//        deform *= sqr((xI-xoff)/norm);
//        this->operator[](ptI) = (amplitude + deform) * sin(omega*t.value());
//        /////////////////////////////////////////////////////////////////////////////////////////////

        // TODO: account for origin not at (0 0 0)

        // get displacement from pre-calculated shape function
        v = vector::zero;
        a = vector::zero;
        for( int inode=0; inode<BD::nnodes; ++inode )
        {
            for( int i=0; i<3; ++i )
            {
                v.component(i) += 
                    BD::h()[ptI*BD::nnodes+inode] * disp[inode].component(i);
                a.component(i) += 
                    BD::h()[ptI*BD::nnodes+inode] * adisp[inode].component(i);
            }
        }


        this->operator[](ptI) = vector::zero;
        // apply rotation
        //for( int j=0; j<3; ++j ) {
        //    for( int i=0; i<3; ++i )
        //    {
        //        this->operator[](ptI).component(j) += R[3*j+i] * v.component(i);
        //    }
        //}

        // TODO: general rotations, retrieve rotation matrix
        // x-rotation
        ang = a.component(0);
        tmp[0] = 0.0;
        tmp[1] = v.component(1)*Foam::cos(ang) - v.component(2)*Foam::sin(ang);
        tmp[2] = v.component(1)*Foam::sin(ang) + v.component(2)*Foam::cos(ang);
        // z-rotation
        ang = a.component(2);
        this->operator[](ptI).component(0) = tmp[0]*Foam::cos(ang) - tmp[1]*Foam::sin(ang);
        this->operator[](ptI).component(1) = tmp[0]*Foam::sin(ang) + tmp[1]*Foam::cos(ang);
        this->operator[](ptI).component(2) = tmp[2];

        if(BD::twoD) this->operator[](ptI).component(BD::bladeDir) = 0.0;

    }

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void beamDynInterfacePointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    beamDynInterfacePointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
