#ifndef DAMAGE_MODEL_H
#define DAMAGE_MODEL_H

class DamageModel
{
	public:
		
		DamageModel();
		
		virtual double getEnvelope(double k_Ang, double dose) = 0;
		virtual int getParameterCount() = 0;
};

#endif
