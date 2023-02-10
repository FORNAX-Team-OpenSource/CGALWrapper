#pragma once

#include "JSONBase.h"
#include <list>
#include <vector>
#include <string>
#include "Point.h"

namespace JSONModels
{
	class Polyline : public JSONBase
	{
	public:
		Polyline();
		virtual ~Polyline();

		virtual bool Deserialize(const rapidjson::Value& obj);
		virtual bool Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const;

		// Getters/Setters.
		std::vector<Point>& PointList() { return _points; }

	private:
		double _area;
		std::vector<Point> _points;

	};
}

