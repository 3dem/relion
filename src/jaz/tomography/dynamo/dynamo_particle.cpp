#include "dynamo_particle.h"
#include <src/jaz/math/Euler_angles_dynamo.h>
#include <sstream>

#define ONE_BASED_COORDS 1
#define OFF_BY_HALF_COORDS 1

using namespace gravis;


DynamoParticle::DynamoParticle()
:	tag(0),
	aligned(0),
	averaged(0),
	dx(0),
	dy(0),
	dz(0),
	tdrot(0),
	tilt(0),
	narot(0),
	cc(0),
	cc2(0),
	cpu(0),
	ftype(0),
	ymintilt(0),
	ymaxtilt(0),
	xmintilt(0),
	xmaxtilt(0),
	fs1(0),
	fs2(0),
	tomo(0),
	reg(0),
	classId(0),
	annotation(0),
	x(0),
	y(0),
	z(0),
	dshift(0),
	daxis(0),
	dnarot(0),
	dcc(0),
	otag(0),
	npar(0),
	ref(0),
	sref(0),
	apix(0),
	def(0),
	eig1(0),
	eig2(0)
{
}

DynamoParticle::DynamoParticle(std::string line)
{
	std::istringstream iss(line);
	
	iss >> tag;
	iss >> aligned;
	iss >> averaged;
	iss >> dx;
	iss >> dy;
	iss >> dz;
	iss >> tdrot;
	iss >> tilt;
	iss >> narot;
	iss >> cc;
	iss >> cc2;
	iss >> cpu;
	iss >> ftype;
	iss >> ymintilt;
	iss >> ymaxtilt;
	iss >> xmintilt;
	iss >> xmaxtilt;
	iss >> fs1;
	iss >> fs2;
	iss >> tomo;
	iss >> reg;
	iss >> classId;
	iss >> annotation;
	iss >> x;
	iss >> y;
	iss >> z;
	iss >> dshift;
	iss >> daxis;
	iss >> dnarot;
	iss >> dcc;
	iss >> otag;
	iss >> npar;
	iss >> ref;
	iss >> sref;
	iss >> apix;
	iss >> def;
	iss >> eig1;
	iss >> eig2;
	
	#if ONE_BASED_COORDS
		x -= 1;
		y -= 1;
		z -= 1;
	#endif
		
	#if OFF_BY_HALF_COORDS
		x += 0.5;
		y += 0.5;
		z += 0.5;
	#endif
		
}

std::string DynamoParticle::getFormattedTag() const
{
	std::string numStr = "xxxxx";
	numStr[4] = '0' + tag % 10;
	numStr[3] = '0' + (tag/10) % 10;
	numStr[2] = '0' + (tag/100) % 10;
	numStr[1] = '0' + (tag/1000) % 10;
	numStr[0] = '0' + (tag/10000) % 10;
	
	return "particle_" + numStr;
}

void DynamoParticle::setAlibiOrientation(d4Matrix A)
{
	x = A(0,3);
	y = A(1,3);
	z = A(2,3);
			
	d3Vector angs = EulerDynamo::matrixToAngles(A);
	
	tdrot = -RAD2DEG(angs[2]);
	tilt  = -RAD2DEG(angs[1]);
	narot = -RAD2DEG(angs[0]);
}

void DynamoParticle::setOrientation(d3Vector north, bool transpose)
{
	north.normalize();
	
	d3Vector left = std::abs(north.x) > std::abs(north.y)? 
				north.cross(d3Vector(0,1,0)).normalize() 
			  : north.cross(d3Vector(1,0,0)).normalize() ;
	d3Vector back = left.cross(north).normalize();
	
	d3Matrix A(
		left.x,  left.y,  left.z,
		back.x,  back.y,  back.z,
		north.x, north.y, north.z );
	
	if (transpose) A.transpose();
		
	const double tol = 1e-16;
	
	if (std::abs(A(2,2)-1.0) < tol)
	{
		narot = 0.0;
		tilt = 0.0;
		tdrot = -atan2(A(1,0),A(0,0));
	}
	else if (std::abs(A(2,2)+1.0) < tol)
	{
		narot = 0;
		tilt  = -PI;
		tdrot = -atan2(A(1,0),A(0,0));
	}
	else
	{
		narot  = -atan2(A(2,0),A(2,1));
		tilt   = -acos(A(2,2));
		tdrot  = -atan2(A(0,2),-A(1,2));
	}
	
	tdrot *= 180.0 / PI;
	tilt  *= 180.0 / PI;
	narot *= 180.0 / PI;
}

std::ostream &operator <<(std::ostream &out, const DynamoParticle &c)
{
    out << c.tag << ' ';
    out << c.aligned << ' ';
    out << c.averaged << ' ';
    out << c.dx << ' ';
    out << c.dy << ' ';
    out << c.dz << ' ';
    out << c.tdrot << ' ';
    out << c.tilt << ' ';
    out << c.narot << ' ';
    out << c.cc << ' ';
    out << c.cc2 << ' ';
    out << c.cpu << ' ';
    out << c.ftype << ' ';
    out << c.ymintilt << ' ';
    out << c.ymaxtilt << ' ';
    out << c.xmintilt << ' ';
    out << c.xmaxtilt << ' ';
    out << c.fs1 << ' ';
    out << c.fs2 << ' ';
    out << c.tomo << ' ';
    out << c.reg << ' ';
    out << c.classId << ' ';
    out << c.annotation << ' ';
	
	double x = c.x;
	double y = c.y;
	double z = c.z;
	
	#if ONE_BASED_COORDS
		x += 1;
		y += 1;
		z += 1;
	#endif
		
	#if OFF_BY_HALF_COORDS
		x -= 0.5;
		y -= 0.5;
		z -= 0.5;
	#endif
	
	out << x << ' ';
	out << y << ' ';
	out << z << ' ';
		
    out << c.dshift << ' ';
    out << c.daxis << ' ';
    out << c.dnarot << ' ';
    out << c.dcc << ' ';
    out << c.otag << ' ';
    out << c.npar << ' ';
    out << c.ref << ' ';
    out << c.sref << ' ';
    out << c.apix << ' ';
    out << c.def << ' ';
    out << c.eig1 << ' ';
    out << c.eig2 << ' ';

    return out;
}
