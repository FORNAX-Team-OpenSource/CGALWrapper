#pragma once

#include "JSONBase.h"
#include <list>
#include <vector>
#include <string>
#include "Polyline.h"

namespace JSONModels
{
	class Polylines : public JSONBase
	{
	public:
		Polylines();
		virtual ~Polylines();

		virtual bool Deserialize(const rapidjson::Value& obj);
		virtual bool Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const;

		// Getters/Setters.
		std::vector<Polyline>& PolylineList() { return _polylines; }

	private:
		std::vector<Polyline> _polylines;

	};
}

