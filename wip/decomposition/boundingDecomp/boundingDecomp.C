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

#include "boundingDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"
#include "globalIndex.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundingDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        boundingDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundingDecomp::boundingDecomp(const dictionary& decompositionDict)
:
    geomDecomp(decompositionDict, typeName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::boundingDecomp::decompose
(
    const polyMesh& polyMesh,
    const pointField& points
)
{
    labelList finalDecomp(points.size());

    if (Pstream::parRun())
    {
        Info << "Foam::labelList Foam::boundingDecomp::decompose : not implemented for parallel" << endl;
        return finalDecomp;
    }

    forAll(finalDecomp, pointI)
    {
        finalDecomp[pointI] = 0;
    }

    label n;
    label target = points.size() / 2;
    const labelList& faceNeighbour = polyMesh.faceNeighbour();
    const labelList& faceOwner = polyMesh.faceOwner();
    for(n = faceNeighbour.size(); n < faceOwner.size(); ++n)
    {
        //Info << "Setting " << faceOwner[n] << " = 1 [ max=" << finalDecomp.size() << " ]" << endl;
        finalDecomp[faceOwner[n]] = 1;
        --target;
    }

    label iter = 1;
    const labelListList& cellCells = polyMesh.cellCells();
    while(target > 0)
    {
        for(label pI = 0; pI < finalDecomp.size() && target > 0; ++pI)
        {
            if(finalDecomp[pI] == iter)
            {
                const labelList& N = cellCells[pI];
                for(label n = 0; n < N.size(); ++n)
                {
                    if(finalDecomp[N[n]] == 0)
                    {
                        finalDecomp[N[n]] = iter+1;
                        //Info << "Setting " << N[n] << " = " << (iter+1) << " [ max=" << finalDecomp.size() << " ]" << endl;
                        --target;
                    }
                }
            }
        }
        ++iter;
    }

    forAll(finalDecomp, pointI)
    {
        if(finalDecomp[pointI] > 0)
        {
            finalDecomp[pointI] = 0;
        }
        else
        {
            finalDecomp[pointI] = 1;
        }
    }

    //Info << "Number of cells: " << points.size() << endl;
    //Info << "Target: " << target << endl;

    return finalDecomp;
}

Foam::labelList Foam::boundingDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cc,
    const scalarField& cWeights,
    labelList& bCells
)
{
    labelList finalDecomp(cc.size(), 0);

    if (Pstream::parRun())
    {
        Info << "Foam::labelList Foam::boundingDecomp::decompose : not implemented for parallel" << endl;
        return finalDecomp;
    }

    label target = cc.size() / 2;

    Info << "Target cells is " << target << endl;
    forAll(bCells, pointI)
    {
        if(bCells[pointI] == 1)
        {
            finalDecomp[pointI] = 1;
            --target;
        }
    }
    Info << "Remaing cells after boundaries is " << target << endl;



    label iter = 1;
    const labelListList& cellCells = globalCellCells;
    while(target > 0)
    {
        for(label pI = 0; pI < finalDecomp.size() && target > 0; ++pI)
        {
            if(finalDecomp[pI] == iter)
            {
                const labelList& N = cellCells[pI];
                for(label n = 0; n < N.size(); ++n)
                {
                    if(finalDecomp[N[n]] == 0)
                    {
                        finalDecomp[N[n]] = iter+1;
                        --target;
                    }
                }
            }
        }
        ++iter;
    }

    forAll(finalDecomp, pointI)
    {
        if(finalDecomp[pointI] > 0)
        {
            finalDecomp[pointI] = 0;
        }
	else
        {
            finalDecomp[pointI] = 1;
        }
    }

    //Info << "Number of cells: " << points.size() << endl;
    //Info << "Target: " << target << endl;

    return finalDecomp;
}



// ************************************************************************* //
