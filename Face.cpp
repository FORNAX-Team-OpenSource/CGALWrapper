#include <string>
#include "Face.h"

namespace JSONModels
{
	Face::Face()
	{}

	Face::~Face()
	{}

	bool Face::Deserialize(const rapidjson::Value& obj)
	{
		//Area(obj["Area"].GetDouble());
		
		for (rapidjson::Value::ConstValueIterator itr = obj["VerticesWithHoles"].Begin(); itr != obj["VerticesWithHoles"].End(); ++itr)
		{
			Polygon p;
			p.Deserialize(*itr);
			_polygons.push_back(p);
		}

		return true;
	}

	bool Face::Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const
	{
		writer->StartObject();

		writer->String("VerticesWithHoles");
		writer->StartArray();

		for (std::vector<Polygon>::const_iterator it = _polygons.begin(); it != _polygons.end(); it++)
		{
			(*it).Serialize(writer);
		}
		writer->EndArray();

		//writer->String("Area");
		//writer->Double(_area);

		writer->EndObject();

		return true;
	}
}
