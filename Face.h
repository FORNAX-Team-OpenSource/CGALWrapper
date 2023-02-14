#pragma once

#include "JSONBase.h"
#include <list>
#include <vector>
#include <string>
#include "Polygon.h"

namespace JSONModels
{
	class Face : public JSONBase
	{
	public:
		Face();
		virtual ~Face();

		virtual bool Deserialize(const rapidjson::Value& obj);
		virtual bool Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const;

		// Getters/Setters.

		//double Area() const { return _area; }
		//void Area(double area) { _area = area; }

		std::vector<Polygon>& PolygonList() { return _polygons; }
		void PolygonList(std::vector<Polygon> polygons) { _polygons = polygons; };

	private:
		//double _area;
		std::vector<Polygon> _polygons;
		
	};
}

