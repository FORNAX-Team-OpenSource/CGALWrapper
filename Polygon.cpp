#include "Polygon.h"

namespace JSONModels
{
	bool Polygon::Deserialize(const rapidjson::Value& obj)
	{
		
		/*rapidjson::Document doc;
		if (!InitDocument(s, doc))
			return false;
			

		if (!doc.IsArray())
			return false;
			*/

		for (rapidjson::Value::ConstValueIterator itr = obj.Begin(); itr != obj.End(); ++itr)
		{
			Point p;
			p.Deserialize(*itr);
			_points.push_back(p);
		}

		return true;
	}

	bool Polygon::Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const
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