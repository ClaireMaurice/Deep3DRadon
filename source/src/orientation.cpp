#include <iomanip>
#include "orientation.h"



/*===========================================================*/
/*                  CLASS MILLER INDICES                     */
/*===========================================================*/
Miller::Miller(int h, int k, int l, int u, int v, int w) : HKL(h,k,l),UVW(u,v,w) {}
Miller::Miller(Eigen::Vector3i hkl, Eigen::Vector3i uvw) : HKL(hkl),UVW(uvw) {}
Miller::Miller(const Miller &mil) : HKL(mil.HKL), UVW(mil.UVW) {} 

Miller &Miller::operator = (const Miller &mil)
{
	HKL = mil.HKL;
	UVW = mil.UVW;
	return (*this);
}

Eigen::Matrix3d Miller::toRotationMatrix()
{
    Eigen::Vector3d hkl;
	hkl[0]=HKL[0] ; hkl[1]=HKL[1] ; hkl[2]=HKL[2];
    hkl.normalize();

    Eigen::Vector3d uvw;
	uvw[0]=UVW[0] ; uvw[1]=UVW[1] ; uvw[2]=UVW[2];
    uvw.normalize();

    Eigen::Vector3d TD = hkl.cross(uvw);

    Eigen::Matrix3d g;
    g.row(0) = uvw;
    g.row(1) = TD;
    g.row(2) = hkl;

    //Eigen::Matrix3d gtrans = g.transpose();

    return g.transpose();
}

Euler Miller::toEuler()
{
	Euler e;

	e.fromRotationMatrix(this->toRotationMatrix());

	return e;
}


std::ostream &operator << (std::ostream &s, const Miller &m) {
    Eigen::Vector3i u,v;
	u = m.Plan();
	v = m.Dir();
	return(s << '(' << u[0] << u[1] << u[2]  << ")[" << v[0] << v[1] << v[2] << ']');
}

/*===========================================================*/
/*                    CLASS EULER ANGLES                     */
/*===========================================================*/
Euler::Euler(double phi1, double PHI, double phi2) : eul(DtoR*phi1, DtoR*PHI, DtoR*phi2) {}
Euler::Euler(const Euler &e) : eul(e.eul) {}
Euler::Euler(const Eigen::Vector3d &vec) : eul(DtoR*vec(0),DtoR*vec(1),DtoR*vec(2)) {}
Euler::Euler(const Eigen::Vector3f &vec) : eul(DtoR*vec(0),DtoR*vec(1),DtoR*vec(2)) {}

Eigen::Matrix3d Euler::toRotationMatrix()
{
    double c1 = cos(eul(0)); double s1 = sin(eul(0));
    double cc = cos(eul(1)) ; double ss = sin(eul(1));
    double c2 = cos(eul(2)) ; double s2 = sin(eul(2));

    Eigen::Matrix3d g;

    g(0,0) = c1*c2 - s1*s2*cc; 	g(0,1) = s1*c2 + c1*s2*cc; 	g(0,2) = s2*ss;
    g(1,0) = -c1*s2 - s1*c2*cc; 	g(1,1) = -s1*s2 + c1*c2*cc; 	g(1,2) = c2*ss;
    g(2,0) = s1*ss; 				g(2,1) = -c1*ss; 				g(2,2) = cc;

	return g;
}

Euler &Euler::fromRotationMatrix(const Eigen::Matrix3d &g)
{
    double f1,F,f2,val;
    val = (fabs(g(2,2))<1.0 ? g(2,2) : g(2,2)>0 ? 1.0 : -1.0);
    F = acos(val);
    if(fabs(g(2,2) - 1.0) > 1.e-6)
	{
        f1 = atan2(g(2,0),-g(2,1));
        f2 = atan2(g(0,2), g(1,2));
	} 
	else 
	{
        f1 = atan2(g(0,1),g(0,0));
		f2 = 0;
	}
    eul << f1,F,f2;
	return (*this);
}

void Euler::asymmetric() {
    // Put PHI in [0-2PI]
    eul[1] -= ((int) (eul[1] / (2 * PI))) * (2 * PI);
    if(eul[1]<0) eul[1]+=2*PI;

    // if F>PI : apply reflection
    if(eul[1]>PI) {
        eul[1] = 2*PI - eul[1];
        eul[0] += PI;
        eul[2] += PI;
    }
    if(eul[1]>1.e-6 && fabs(eul[1]-PI)>1.e-6) {
        /* ranging phi2 within [0,2PI] */
        eul[2] -= ((int) (eul[2] / (2 * PI))) * (2 * PI);
        if (eul[2] < 0)
          eul[2] += 2*PI;
    } else if( eul[1]>PI )
    { eul[0] += eul[2]; eul[2]=0;}
    else {eul[0] -= eul[2]; eul[2]=0;}

    eul[0] -= ((int) (eul[0] / (2*PI))) * (2*PI);
    if (eul[0] < 0)
      eul[0] += (2 * PI);
}

Euler Euler::get_fromRotationMatrix(const Eigen::Matrix3d &g) {
    Euler e;
    e.fromRotationMatrix(g);
    return e;
}

Euler &Euler::smallestVariant(int spaceGroup) {
    Quat qm(*this);

    if(spaceGroup >= 195) {
        //qDebug() << this->phi1() << " " << this->phi() << " " << this->phi2();
        double romax=0;
        for(int iSym=0 ; iSym<24 ; ++iSym)
        {
            Quat qi;
            qi = qm * SymQuatCubic[iSym];
            if(qi.ro()>romax) {
                *this = qi.toEuler();
                romax=qi.ro();
            }

        }
    }


    return *this;
}

std::ostream &operator << (std::ostream &s, const Euler &e) {
	return(s << '(' << e.phi1() << ", " << e.phi() << ", " << e.phi2()  << ')' );
}


/*===========================================================*/
/*                    CLASS QUATERNION                       */
/*===========================================================*/
Quat::Quat(double ro, double lambda, double mu, double nu) : scal(ro),vec(lambda,mu,nu) {}

Quat::Quat(double ro, Eigen::Vector3d v, bool AngleAxis)
{
	if(AngleAxis)
	{
		scal = cos(DtoR * ro / 2);
        vec = v.normalized() * sin(DtoR * ro / 2);
	}
	else
	{	
		scal = ro;
		vec = v;
	}
}

Quat::Quat(const Quat &q) : scal(q.scal),vec(q.vec) {}

Quat::Quat(const Euler &e) {
    scal = cos(e.phi()*DtoR/2)*cos( (e.phi1()+e.phi2())*DtoR/2);
    vec[0] = sin(e.phi()*DtoR/2)*cos( (e.phi1()-e.phi2())*DtoR/2);
    vec[1] = sin(e.phi()*DtoR/2)*sin( (e.phi1()-e.phi2())*DtoR/2);
    vec[2] = cos(e.phi()*DtoR/2)*sin( (e.phi1()+e.phi2())*DtoR/2);
}

Quat::Quat(const Rodrigues &R) {
    Eigen::Vector3d r(R.R1(),R.R2(),R.R3());
    double tanOmega = r.norm();
    scal=cos(atan(tanOmega));
    if(tanOmega>1.e-6)
        vec=sin(atan(tanOmega))*r.normalized();
    else
        vec<<0,0,0;
}

Quat::Quat(const AxisAnglePair &aap) {
    double angle(aap.Angle()*DtoR/2);
    scal = cos(angle);
    vec = sin(angle)*aap.Axis();
}

Quat::Quat(const double theta, const int axisindex) : 
      scal(cos(DtoR * theta / 2)),
      vec((axisindex==0 ? sin(DtoR * theta / 2) : 0.0),
          (axisindex==1 ? sin(DtoR * theta / 2) : 0.0),
          (axisindex==2 ? sin(DtoR * theta / 2) : 0.0))
{
}


double &Quat::operator[](int index) 
{
	if(index==0) return (scal);
	return (vec[index-1]);
	//return(0);
}
const double &Quat::operator[](int index) const
{
	if(index==0) return (scal);
    return (this->vec.coeff(index-1));
	//return(0);
}

Quat &Quat::operator = (const Quat &q) 
{
	scal = q.ro();
	vec = q.v();
	return(*this);
}

Quat &Quat::operator += (const Quat &q) 
{
	scal += q.ro();
	vec += q.v();
	return(*this);
}

Quat &Quat::operator -= (const Quat &q)
{
	scal -= q.ro();
	vec -= q.v();
	return(*this);
}

Quat &Quat::operator *= (double alpha)
{
	scal *= alpha;
	vec *= alpha;
	return(*this);
}
Quat &Quat::operator *= (const Quat &q)			// attention : non commutatif
{
    double newscal = scal*q.ro() - vec.dot(q.v());
    Eigen::Vector3d newvec = scal*q.v() + q.ro()*vec + vec.cross(q.v());

	if(newscal < 0) {newscal *= -1 ; newvec *= -1;}
	scal = newscal;
	vec = newvec;
	
	return(*this);
}
Quat &Quat::operator /= (double alpha)
{
	scal /= alpha;
	vec /= alpha;
	return(*this);
}

// comparaison
bool Quat::operator == (const Quat &q) const
{
	return (scal==q.ro() && vec==q.v());
}
bool Quat::operator != (const Quat &q) const
{
	return( scal != q.ro() || vec != q.v());
}

// operateurs arithmetiques
Quat Quat::operator + (const Quat &q) const
{
	Quat result;

	result.scal = scal + q.scal;
	result.vec = vec + q.vec;

	return result;
}

Quat Quat::operator - (const Quat &q) const
{
	Quat result;

	result.scal = scal - q.scal;
	result.vec = vec - q.vec;

	return result;
}


Quat Quat::operator - () const 
{
	Quat result;

	result.scal = -scal;
	result.vec = -vec;

	return result;
}

Quat Quat::operator * (const Quat &q) const
{
	Quat result;

    result.scal = scal*q.ro() - vec.dot(q.v());
    result.vec = scal*q.v() + q.ro()*vec + vec.cross(q.v());

	if(result.scal < 0) result = -result;

	return result;
}

Quat &Quat::Normalise()
{
	double norm;
    norm = sqrt(scal*scal + vec.squaredNorm());
	(*this) /= norm;
	return (*this);
}

Quat &Quat::Conjugate()
{
	vec = -vec;
	return (*this);
}

Quat &Quat::Inverse()
{
	vec = -vec;					// uniquement pour quaternions unitaires
	return (*this);				// TO DO : fiabiliser pour quaternions quelconques
}

Quat Quat::getInverse()
{
    return Quat(*this).Inverse();
}

Eigen::Matrix3d Quat::toRotationMatrix()
{
    Eigen::Matrix3d g;
	double ro( (*this).ro()), lambda( (*this).lambda()), mu( (*this).mu()), nu( (*this).nu());
    double ro2(ro*ro),lambda2(lambda*lambda),mu2(mu*mu),nu2(nu*nu);

    g(0,0) = ro2 + lambda2 - mu2 - nu2;
    g(1,1) = ro2 - lambda2 + mu2 - nu2;
    g(2,2) = ro2 - lambda2 - mu2 + nu2;

    //double check1=g(0,0)+g(1,1)+g(2,2)+1;
    //double check2=4*ro2;

    g(0,1) = g(1,0) = lambda*mu;
    g(0,1) += ro*nu ; g(1,0) -= ro*nu;
    g(0,1) *= 2; g(1,0) *= 2;

    g(0,2) = g(2,0) = lambda*nu;
    g(0,2) -= ro*mu; ; g(2,0) += ro*mu;
    g(0,2) *= 2; g(2,0) *= 2;
	
    g(1,2) = g(2,1) = mu*nu;
    g(1,2) += ro*lambda;	g(2,1) -= ro*lambda;
    g(1,2) *= 2; g(2,1) *= 2;

	return g;
}

Quat &Quat::fromRotationMatrix(const Eigen::Matrix3d &g, bool sym, int spacegroup)
{
    Quat q;
    q.scal = 0.5*sqrt(g(0,0)+g(1,1)+g(2,2) + 1);
    if(q.scal > 1.e-6)
	{
        q.vec[0] = (g(1,2) - g(2,1))/q.scal/4;
        q.vec[1] = (g(2,0) - g(0,2))/q.scal/4;
        q.vec[2] = (g(0,1) - g(1,0))/q.scal/4;
	}
	else		// ro=0 : half-turn
	{
        q.vec[0] = sqrt( (g(0,0) + 1)/2 );
        q.vec[1] = sqrt( (g(1,1) + 1)/2 );
        q.vec[2] = sqrt( (g(2,2) + 1)/2 );

		// attribution des signes : par convention le terme le plus grand (m) est >0
		// les autres ont le m\EAme signe que g[i][m]
        int m(-1);
        for(int i=0 ; i<3 ; ++i)
		{
            if(q.vec[i]>=q.vec[(i+1)%3] && q.vec[i]>=q.vec[(i+2)%3])
			{
				m=i; break;
			}
		}
        for(int i=0 ; i<3 ; ++i)
		{
            if(i!=m) q.vec[i] *= g(i,m)/fabs(g(i,m));
		}
	}

    q.Normalise();

    if(sym) {
        // on cherche le quaternion correspondant a la rotation equivalente (au sens de la symetrie cubique)
        // d'angle minimum (partie scalaire maximum)
        Quat SymQuat[24],qi;
        double rac2 = 1/sqrt(2.0);

        SymQuat[0].scal = 1; SymQuat[0].vec << 0,0,0;

        SymQuat[1].scal = rac2; SymQuat[1].vec << rac2,0,0;
        SymQuat[2].scal = 0; SymQuat[2].vec << 1,0,0;
        SymQuat[3].scal = -rac2; SymQuat[3].vec << rac2,0,0;

        SymQuat[4].scal = rac2; SymQuat[4].vec << 0,rac2,0;
        SymQuat[5].scal = 0; SymQuat[5].vec << 0,1,0;
        SymQuat[6].scal = -rac2; SymQuat[6].vec << 0,rac2,0;

        SymQuat[7].scal = rac2; SymQuat[7].vec << 0,0,rac2;
        SymQuat[8].scal = 0; SymQuat[8].vec << 0,0,1;
        SymQuat[9].scal = -rac2; SymQuat[9].vec << 0,0,rac2;

        SymQuat[10].scal = 0; SymQuat[10].vec << rac2,rac2,0;
        SymQuat[11].scal = 0; SymQuat[11].vec << -rac2,rac2,0;

        SymQuat[12].scal = 0; SymQuat[12].vec << 0,rac2,rac2;
        SymQuat[13].scal = 0; SymQuat[13].vec << 0,-rac2,rac2;

        SymQuat[14].scal = 0; SymQuat[14].vec << rac2,0,rac2;
        SymQuat[15].scal = 0; SymQuat[15].vec << -rac2,0,rac2;

        SymQuat[16].scal = 0.5; SymQuat[16].vec << 0.5,0.5,0.5;
        SymQuat[17].scal = 0.5; SymQuat[17].vec << -0.5,0.5,0.5;
        SymQuat[18].scal = 0.5; SymQuat[18].vec << 0.5,-0.5,0.5;
        SymQuat[19].scal = 0.5; SymQuat[19].vec << 0.5,0.5,-0.5;

        SymQuat[20].scal = -0.5; SymQuat[20].vec << 0.5,0.5,0.5;
        SymQuat[21].scal = -0.5; SymQuat[21].vec << -0.5,0.5,0.5;
        SymQuat[22].scal = -0.5; SymQuat[22].vec << 0.5,-0.5,0.5;
        SymQuat[23].scal = -0.5; SymQuat[23].vec << 0.5,0.5,-0.5;

        double SCAL=0;
        for(int i=0 ; i<24 ; ++i) {
            SymQuat[i].Normalise();
            qi = q*SymQuat[i];
            Eigen::Matrix3d Ri;
            Ri = qi.toRotationMatrix();
            if(qi.scal > SCAL) {
                SCAL = qi.scal, *this=qi;
            }
        }

    } else
        *this=q;

	return (*this);
}

Quat & Quat::MinimalMisorientationQuat(const Quat &q1, const Quat &q2,
                                       int spacegroup, bool active)
{
    Quat invq1(q1),qm_crystal,qm_sample;
    invq1.Inverse();

    if(spacegroup >= 195)  {  // symetrie cubique

//       Quat SymQuat[24],qi;
//       double rac2 = 1/sqrt(2.0);

//       SymQuat[0].scal = 1; SymQuat[0].vec << 0,0,0;
//       SymQuat[1].scal = rac2; SymQuat[1].vec << rac2,0,0;
//       SymQuat[2].scal = 0; SymQuat[2].vec << 1,0,0;
//       SymQuat[3].scal = -rac2; SymQuat[3].vec << rac2,0,0;

//       SymQuat[4].scal = rac2; SymQuat[4].vec << 0,rac2,0;
//       SymQuat[5].scal = 0; SymQuat[5].vec << 0,1,0;
//       SymQuat[6].scal = -rac2; SymQuat[6].vec << 0,rac2,0;

//       SymQuat[7].scal = rac2; SymQuat[7].vec << 0,0,rac2;
//       SymQuat[8].scal = 0; SymQuat[8].vec << 0,0,1;
//       SymQuat[9].scal = -rac2; SymQuat[9].vec << 0,0,rac2;

//       SymQuat[10].scal = 0; SymQuat[10].vec << rac2,rac2,0;
//       SymQuat[11].scal = 0; SymQuat[11].vec << -rac2,rac2,0;

//       SymQuat[12].scal = 0; SymQuat[12].vec << 0,rac2,rac2;
//       SymQuat[13].scal = 0; SymQuat[13].vec << 0,-rac2,rac2;

//       SymQuat[14].scal = 0; SymQuat[14].vec << rac2,0,rac2;
//       SymQuat[15].scal = 0; SymQuat[15].vec << -rac2,0,rac2;

//       SymQuat[16].scal = 0.5; SymQuat[16].vec << 0.5,0.5,0.5;
//       SymQuat[17].scal = 0.5; SymQuat[17].vec << -0.5,0.5,0.5;
//       SymQuat[18].scal = 0.5; SymQuat[18].vec << 0.5,-0.5,0.5;
//       SymQuat[19].scal = 0.5; SymQuat[19].vec << 0.5,0.5,-0.5;

//       SymQuat[20].scal = -0.5; SymQuat[20].vec << 0.5,0.5,0.5;
//       SymQuat[21].scal = -0.5; SymQuat[21].vec << -0.5,0.5,0.5;
//       SymQuat[22].scal = -0.5; SymQuat[22].vec << 0.5,-0.5,0.5;
//       SymQuat[23].scal = -0.5; SymQuat[23].vec << 0.5,0.5,-0.5;

        double SCAL=0;

        for(int iSym=0 ; iSym<24 ; ++iSym)
        {
            //qm_crystal = invq1 * q2*SymQuatCubic[iSym];  // in crystal1 coordinate system
            qm_sample = q2*SymQuatCubic[iSym]*invq1;    // in sample coordinate system
            if(qm_sample.scal > 1.0)  // treat numerical errors when identical orientations
            {
                SCAL=1.0;
                qm_sample.scal=1.0;
                *this = qm_sample;
            }
            if(qm_sample.scal>SCAL) {
                SCAL = qm_sample.scal;
                *this = qm_sample;
            }
            if(SCAL > cos(30*PI/360)) break;
        }
    } else {
        if(spacegroup>=168) {       // symetrie hexagonale - 12 elements de sym
            double SCAL=0;

            for(int iSym=0 ; iSym<12 ; ++iSym)
            {
               qm_sample = q2*SymQuatHexagonal[iSym]*invq1;
               if(qm_sample.scal > 1.0)
               {
                  SCAL=1.0;
               }
               if(qm_sample.scal>SCAL) {
                  SCAL = qm_sample.scal;
                  *this = qm_sample;
               }
            }
        }
    }
    return *this;
}

Quat &Quat::sortDecreasingOrder()
{
    Eigen::Vector4d v(fabs(ro()),fabs(lambda()),fabs(mu()),fabs(nu()));

    int imin,imax;
    v.minCoeff(&imin);

    double swap;
    swap = v[3];v[3]=v[imin];v[imin]=swap;
    v.maxCoeff(&imax);
    swap = v[0];v[0]=v[imax];v[imax]=swap;
    if(v[1]<v[2]) {
        swap = v[1];v[1]=v[2];v[2]=swap;
    }
    scal = v[0];
    vec = v.block(1,0,3,1);
    return (*this);				// TO DO : fiabiliser pour quaternions quelconques
}

double Quat::desorientationAngle(const Quat &q2, int spaceGroup) {

    if(spaceGroup >= 195) {
        Quat invq1(*this),qm;

        invq1.Inverse();
        qm = invq1 * q2;

        double angle=0;
        qm=qm.sortDecreasingOrder();
        Eigen::Vector3d variants(qm[0],(qm[0]+qm[1])/sqrt(2.0),(qm[0]+qm[1]+qm[2]+qm[3])/2.0);
        angle = variants.maxCoeff();

        return (angle>=1.0 ? 0 : RtoD*2*acos(angle));
    } else {
        Quat qm=MinimalMisorientationQuat(*this,q2,spaceGroup);
        return (fabs(qm.ro())>=1.0 ? 0 : RtoD*2*acos(qm.ro()));
    }

    return 0.0;

}


Euler Quat::toEuler()
{
	Euler e;

	// calcul de PHI
    e.set_phi(2*atan2(sqrt(vec[0]*vec[0] + vec[1]*vec[1]),
                      sqrt(scal*scal + vec[2]*vec[2])));


	if(e.phi()!=0 && e.phi() != 180)
	{
		e.set_phi1(atan2(vec[2],scal) + atan2(vec[1],vec[0]));
		e.set_phi2(atan2(vec[2],scal) - atan2(vec[1],vec[0]));
	} 
	else
	{
		if(e.phi()==0)
		{
            e.set_phi1(2*atan2(vec[2],scal));
			e.set_phi2(0);
		}
		if(e.phi()==180)
		{
            e.set_phi1(2*atan2(vec[1],vec[0]));
			e.set_phi2(0);
		}
	}
	return e;
}

Quat &Quat::fromEuler(const Euler &e)
{
    double fi(e.phi()*DtoR/2);
    double fi1pfi2((e.phi1()+e.phi2())*DtoR/2);
    double fi1mfi2((e.phi1()-e.phi2())*DtoR/2);

    scal = cos(fi)*cos(fi1pfi2);
    vec[0] = sin(fi)*cos(fi1mfi2);
    vec[1] = sin(fi)*sin(fi1mfi2);
    vec[2] = cos(fi)*sin(fi1pfi2);

	return (*this);
}

Quat &Quat::fromAxisAnglePair(const AxisAnglePair &aap)
{
	double angle(aap.Angle()*DtoR/2);
	scal = cos(angle);
	vec = sin(angle)*aap.Axis();

	return (*this);
}

Quat &Quat::fromAngleAxis(const double theta, const int axisindex)
{
   double angle=DtoR*theta/2.0;
   scal=cos(angle);
   vec << (axisindex==0 ? sin(angle) : 0.0),
          (axisindex==1 ? sin(angle) : 0.0),
          (axisindex==2 ? sin(angle) : 0.0);
	return (*this);
}

std::ostream &operator << (std::ostream &s, const Quat &q) {
	return(s << '[' << q.ro() << ", " << q.v() <<  ']');
}


/*===========================================================*/
/*                  CLASS RODRIGUES VECTOR                   */
/*===========================================================*/
Rodrigues::Rodrigues(const double R1, const double R2, const double R3) : R(R1,R2,R3) {}
Rodrigues::Rodrigues(const Eigen::Vector3d &v) : R(v) {}
Rodrigues::Rodrigues(const Rodrigues &A) : R(A.R) {}
Rodrigues::Rodrigues(const AxisAnglePair AAP)
{
    R = tan(0.5*DtoR*AAP.Angle()) * AAP.Axis();
}

int permutation(int i,int j,int k)
{
    if(i==0 && j==1 && k==2) return 1;
    if(i==1 && j==2 && k==0) return 1;
    if(i==2 && j==0 && k==1) return 1;

    if(i==0 && j==2 && k==1) return -1;
    if(i==1 && j==0 && k==2) return -1;
    if(i==2 && j==1 && k==0) return -1;

    return 0;
}

Rodrigues &Rodrigues::operator = (const Rodrigues &Rod)
{
    R = Rod.R;
    return(*this);
}

Rodrigues &Rodrigues::fromEigenVector(const Eigen::Vector3d &v) {
    this->R = v;
    return(*this);
}

Rodrigues &Rodrigues::fromQuaternion(const Quat &q)
{
    R << q.lambda()/q.ro(),q.mu()/q.ro(), q.nu()/q.ro();
    return(*this);
}

Quat &Quat::fromRodrigues(const Rodrigues &R) {
    double normedeR = R.vec().norm();
    double thetasur2 = atan(R.vec().norm());
    (*this).scal = cos(thetasur2);
    (*this).vec = R.vec()/sin(thetasur2);
    return(*this);
}

Eigen::Matrix3d Rodrigues::toRotationMatrix() {
    Eigen::Matrix3d ROT;

    double ro_lxro_l = R.squaredNorm();
    for(int i=0 ; i<3 ; ++i)
        for(int j=0 ; j<3 ; ++j) {
            ROT(i,j) = ( (1-ro_lxro_l)*(i==j ? 1 : 0) + 2*R(i)*R(j) );
            for(int k=0 ; k<3 ; ++k)
                ROT(i,j) += 2*permutation(i,j,k)*R(k);
            ROT(i,j) /= (1+ro_lxro_l);
        }

    return ROT;
}



/*===========================================================*/
/*                 CLASS AXIS/ANGLE PAIR                     */
/*===========================================================*/
AxisAnglePair::AxisAnglePair(double r1, double r2, double r3, double angle) : r(r1,r2,r3),theta(DtoR*angle)
{
    r.normalize();
}

AxisAnglePair::AxisAnglePair(const Eigen::Vector3d &axis, double angle) : r(axis),theta(DtoR*angle)
{
    r.normalize();
}

AxisAnglePair::AxisAnglePair(const AxisAnglePair &AAP) : r(AAP.r),theta(AAP.theta) {}

AxisAnglePair &AxisAnglePair::operator = (const AxisAnglePair &AAP) 
{
	theta = AAP.theta;
	r = AAP.r;
	return(*this);
}

AxisAnglePair &AxisAnglePair::fromRotationMatrix(const Eigen::Matrix3d g)
{
    double trace=(g.trace()-1.0)/2.0;
    theta = acos( trace<=1.0 ? trace : 1.0 );
	
	if(theta >= 1.e-6 && fabs(theta - PI) >= 1.e-6)
	{
        r[0] = (g(1,2)-g(2,1));
        r[1] = (g(2,0)-g(0,2));
        r[2] = (g(0,1)-g(1,0));
		double val=2*sin(theta);
		r /= val;
	}
	else
	{
		if(theta < 1.e-6) {
			r[0]=1 ; r[2]=r[1]=0 ;
		}
		else
		{
            r[0] = sqrt( (g(0,0)+1)/2 );
            r[1] = sqrt( (g(1,1)+1)/2 );
            r[2] = sqrt( (g(2,2)+1)/2 );

			// attribution des signes : par convention le terme le plus grand (m) est >0
			// les autres ont le m\EAme signe que g[i][m]
            int m(-1);
            for(int i=0 ; i<3 ; ++i)
			{
				if(r[i]>=r[(i+1)%3] && r[i]>=r[(i+2)%3]) 
				{
					m=i; break;
				}
			}
            for(int i=0 ; i<3 ; ++i)
			{
                if(i!=m) r[i] *= g(i,m)/fabs(g(i,m));
			}

		}
	}
    r.normalize();
	return (*this);
}

Eigen::Matrix3d AxisAnglePair::toRotationMatrix()
{
	#define r1 r[0]
	#define r2 r[1]
	#define r3 r[2]

    Eigen::Matrix3d g;
	double c(cos(theta)),s(sin(theta));

    g(0,0) = r1*r1*(1-c) + c;
    g(0,1) = r1*r2*(1-c) + r3*s;
    g(0,2) = r1*r3*(1-c) - r2*s;
    g(1,0) = r2*r1*(1-c) - r3*s;
    g(1,1) = r2*r2*(1-c) + c;
    g(1,2) = r2*r3*(1-c) + r1*s;
    g(2,0) = r3*r1*(1-c) + r2*s;
    g(2,1) = r3*r2*(1-c) - r1*s;
    g(2,2) = r3*r3*(1-c) + c;

	return g;
}

