#pragma once

#include "JSONBase.h"
#include <list>
#include <string>

namespace JSONModels
{
	class Point : public JSONBase
	{
	public:
		Point();
		virtual ~Point();

		virtual bool Deserialize(const rapidjson::Value& obj);
		virtual bool Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const;

		// Getters/Setters.
		float X() const { return _x; }
		void X(double x) { _x = x; }

		float Y() const { return _y; }
		void Y(double y) { _y = y; }

		float Z() const { return _z; }
		void Z(double z) { _z = z; }
	private:
		double _x;
		double _y;
		double _z;
	};
}

