#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Common.h"

namespace crpropa {

PeriodicMagneticField::PeriodicMagneticField(ref_ptr<MagneticField> field,
		const Vector3d &extends) :
		field(field), extends(extends), origin(0, 0, 0), reflective(false) {

}

PeriodicMagneticField::PeriodicMagneticField(ref_ptr<MagneticField> field,
		const Vector3d &extends, const Vector3d &origin, bool reflective) :
		field(field), extends(extends), origin(origin), reflective(reflective) {

}

Vector3d &PeriodicMagneticField::getOrigin() {
	return origin;
}

void PeriodicMagneticField::setOrigin(const Vector3d &origin) {
	this->origin = origin;
}

Vector3d &PeriodicMagneticField::getExtends() {
	return extends;
}

void PeriodicMagneticField::setExtends(const Vector3d &origin) {
	this->extends = extends;
}

bool PeriodicMagneticField::isReflective() {
	return reflective;
}

void PeriodicMagneticField::setReflective(bool reflective) {
	this->reflective = reflective;
}

Vector3d PeriodicMagneticField::getField(const Vector3d &position) const {
	Vector3d n = ((position - origin) / extends).floor();
	Vector3d p = position - origin - n * extends;

	if (reflective) {
		long mx = (long) ::fabs(n.x) % 2;
		if (mx == 1)
			p.x = extends.x - p.x;
		long my = (long) ::fabs(n.y) % 2;
		if (my == 1)
			p.y = extends.y - p.y;
		long mz = (long) ::fabs(n.z) % 2;
		if (mz == 1)
			p.z = extends.z - p.z;
	}

	return field->getField(p);
}

void MagneticFieldList::addField(ref_ptr<MagneticField> field) {
	fields.push_back(field);
}

Vector3d MagneticFieldList::getField(const Vector3d &position) const {
	Vector3d b;
	for (int i = 0; i < fields.size(); i++)
		b += fields[i]->getField(position);
	return b;
}

MagneticFieldEvolution::MagneticFieldEvolution(ref_ptr<MagneticField> field,
	double m) :
	field(field), m(m) {
}

Vector3d MagneticFieldEvolution::getField(const Vector3d &position,
	double z) const {
	return field->getField(position) * pow(1+z, m);
}

Vector3d MagneticDipoleField::getField(const Vector3d &position) const {
		Vector3d r = (position - origin);
		Vector3d unit_r = r.getUnitVector();
		
		if (r.getR() == 0) { // singularity
			return moment * 2 * mu0 / 3;
		}
		return (unit_r * (unit_r.dot(moment)) * 3 - moment) / pow(r.getR() / radius, 3) * mu0 / (4*M_PI);
}

RadialMagneticField::RadialMagneticField(ref_ptr<MagneticField> field, const Vector3d &origin, const double r0) :
	        field(field) {
    setScaleRadius(r0);
	setOrigin(origin);
}

//rn only sets origin to (0,0,0), but could maybe be expanded to get origin from original magnetic field if applicable (probably more an edge case )
RadialMagneticField::RadialMagneticField(ref_ptr<MagneticField> field, const double r0) :
	        field(field) {
    setScaleRadius(r0);
	setOrigin(Vector3d(0.)); //is this legal?
}

Vector3d RadialMagneticField::getField(const Vector3d &position, double z, double t) const {
	Vector3d r = (position - origin);
	Vector3d B(0,0,0);
	if (field.valid()) {
		B = field->getField(position, z, t);
		B *= pow(1 / (1 + pow_integer<2>(r.getR() / r0)), 2./3.);
	}

	return B;
}

void RadialMagneticField::setField(ref_ptr<MagneticField> field) {
    this->field = field;
}
ref_ptr<MagneticField> RadialMagneticField::getField() const {
    return field;
}

Vector3d RadialMagneticField::getOrigin() {
	return origin;
}
void RadialMagneticField::setOrigin(const Vector3d &origin) {
	this->origin = origin;
}

void RadialMagneticField::setScaleRadius(const double r0) {
	this->r0 = r0;
}
double RadialMagneticField::getScaleRadius() const {
	return r0;
}

#ifdef CRPROPA_HAVE_MUPARSER
RenormalizeMagneticField::RenormalizeMagneticField(ref_ptr<MagneticField> field,
		std::string expression) :
		field(field), expression(expression) {

	p =  new mu::Parser();
	p->DefineVar("B", &Bmag);
	p->DefineConst("tesla", tesla);
	p->DefineConst("gauss", gauss);
	p->DefineConst("muG", muG);
	p->DefineConst("nG", nG);
	p->SetExpr(expression);
}

Vector3d RenormalizeMagneticField::getField(const Vector3d &position) {
	Vector3d B = field->getField(position);
	Bmag = B.getR();
	return B * p->Eval();
}
#endif

} // namespace crpropa
