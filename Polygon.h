#pragma once

#include "Point.h"
#include <list>
#include <vector>
#include <string>

namespace JSONModels
{
	class Polygon : public JSONBase
	{
	public:
		virtual ~Polygon() {};

		std::vector<Point>& PointsList() { return _points; }
		void PointsList(std::vector<Point> points) { _points = points; }

		virtual bool Deserialize(const rapidjson::Value& obj);
		virtual bool Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const;
	private:
		std::vector<Point> _points;
	};
}

