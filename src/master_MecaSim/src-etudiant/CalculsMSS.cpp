/*
 * CalculsMSS.cpp :
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/** \file CalculsMSS.cpp
Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant 
 (methode d 'Euler semi-implicite) : principales fonctions de calculs.
\brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
*/ 

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleMSS.h"
#include "Viewer.h"

using namespace std;





/**
* Calcul des forces appliquees sur les particules du systeme masses-ressorts.
 */
void ObjetSimuleMSS::CalculForceSpring()
{
	/// f = somme_i (ki * (l(i,j)-l_0(i,j)) * uij ) + (nuij * (vi - vj) * uij) + (m*g) + force_ext
	
	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration

    vector<Ressort*> ressorts = _SytemeMasseRessort->GetRessortList();

	for (int i = 0; i < _SytemeMasseRessort->GetNbRessort(); i++) {
	    Particule *pA = ressorts[i]->GetParticuleA();
	    Particule *pB = ressorts[i]->GetParticuleB();

	    float longeur = length(P[pA->GetId()] - P[pB->GetId()]);

	    Vector uijA = (P[pA->GetId()] - P[pB->GetId()]) / longeur;
	    Vector uijB = (P[pB->GetId()] - P[pA->GetId()]) / longeur;

	    Vector elaA = ressorts[i]->GetRaideur() * (longeur - ressorts[i]->GetLrepos()) * uijB;
	    Vector elaB = ressorts[i]->GetRaideur() * (longeur - ressorts[i]->GetLrepos()) * uijA;

	    Vector viscoA = (ressorts[i]->GetAmortissement() * dot(V[pB->GetId()] - V[pA->GetId()], uijB)) * uijB;
	    Vector viscoB = (ressorts[i]->GetAmortissement() * dot(V[pA->GetId()] - V[pB->GetId()], uijA)) * uijA;

	    Force[pA->GetId()] = Force[pA->GetId()] + elaA + viscoA;
        Force[pB->GetId()] = Force[pB->GetId()] + elaB + viscoB;
	}
    
		
}//void


/**
 * Gestion des collisions avec le sol.
 */
void ObjetSimuleMSS::CollisionPlan()
{
    /// Arret de la vitesse quand touche le plan
	vector<Particule*> parts = _SytemeMasseRessort->GetPartList();
	for (int i = 0; i < _SytemeMasseRessort->GetNbParticule(); i++) {
		if(P[i].y <= -10) {
		    P[i].y = -20 - P[i].y;
		    V[i].y = -V[i].y * 0.9;
		}
	}

}// void

/**
 * Gestion de la collision avec un autre objet
 */
void ObjetSimuleMSS::Collision(Transform t) {

}//void

void ObjetSimuleMSS::CollisionSphere(Point p, double rayon, double visco) {
    for (int i = 0; i < P.size(); i++) {
        if( distance(P[i] + Point(0,0,0), p) < rayon) {
            Vector cp(p,P[i] + Point(0,0,0));
            Vector n = normalize(cp);

            P[i] = Vector(Point(0,0,0) ,(n * rayon + p));
            V[i]= visco * (V[i] + 2 * dot(V[i], -n) * n);
        }
    }
}
