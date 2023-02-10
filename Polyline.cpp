#include <string>
#include "Polyline.h"

namespace JSONModels
{
	Polyline::Polyline()
	{}

	Polyline::~Polyline()
	{}

	bool Polyline::Deserialize(const rapidjson::Value& obj)
	{
	
		//for (rapidjson::Value::ConstValueIterator itr = obj["VerticesWithHoles"].Begin(); itr != obj["VerticesWithHoles"].End(); ++itr)
		//{
		//	Point p;
		//	p.Deserialize(*itr);
		//	_points.push_back(p);
		//}

		return true;
	}

	bool Polyline::Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const
	{
		writer->StartArray();
		for (std::vector<Point>::const_iterator it = _points.begin(); it != _points.end(); it++)
		{
			(*it).Serialize(writer);
		}
		writer->EndArray();
		return true;
	}
}
