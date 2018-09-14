// License: Apache 2.0. See LICENSE file in root directory.
// Copyright(c) 2015 Intel Corporation. All Rights Reserved.

#pragma once

namespace rs2 {

	class color_map
	{
	public:
		color_map(std::map<float, float3> map, int steps = 4000) : _map(map)
		{
			initialize(steps);
		}

		color_map(const std::vector<float3>& values, int steps = 4000)
		{
			for (size_t i = 0; i < values.size(); i++)
			{
				_map[(float)i / (values.size() - 1)] = values[i];
			}
			initialize(steps);
		}

		color_map() {}

		inline float3 get(float value) const
		{
			if (_max == _min) return *_data;
			auto t = (value - _min) / (_max - _min);
			t = clamp(t, 0.f, 1.f);
			return _data[(int)(t * (_size - 1))];
		}

		float min_key() const { return _min; }
		float max_key() const { return _max; }

	private:
		inline float3 lerp(const float3& a, const float3& b, float t) const
		{
			return b * t + a * (1 - t);
		}

		float3 calc(float value) const
		{
			if (_map.size() == 0) return{ value, value, value };
			if (_map.find(value) != _map.end()) return _map.at(value);

			if (value < _map.begin()->first)   return _map.begin()->second;
			if (value > _map.rbegin()->first)  return _map.rbegin()->second;

			auto lower = _map.lower_bound(value) == _map.begin() ? _map.begin() : --(_map.lower_bound(value));
			auto upper = _map.upper_bound(value);

			auto t = (value - lower->first) / (upper->first - lower->first);
			auto c1 = lower->second;
			auto c2 = upper->second;
			return lerp(c1, c2, t);
		}

		void initialize(int steps)
		{
			if (_map.size() == 0) return;

			_min = _map.begin()->first;
			_max = _map.rbegin()->first;

			_cache.resize(steps + 1);
			for (int i = 0; i <= steps; i++)
			{
				auto t = (float)i / steps;
				auto x = _min + t*(_max - _min);
				_cache[i] = calc(x);
			}

			_size = _cache.size();
			_data = _cache.data();
		}

		std::map<float, float3> _map;
		std::vector<float3> _cache;
		float _min, _max;
		size_t _size; float3* _data;
	};

	static color_map jet{ {
		{ 0, 0, 255 },
		{ 0, 255, 255 },
		{ 255, 255, 0 },
		{ 255, 0, 0 },
		{ 50, 0, 0 },
		} };

	static color_map classic{ {
		{ 30, 77, 203 },
		{ 25, 60, 192 },
		{ 45, 117, 220 },
		{ 204, 108, 191 },
		{ 196, 57, 178 },
		{ 198, 33, 24 },
		} };

	static color_map grayscale{ {
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		} };

	static color_map inv_grayscale{ {
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		} };

	static color_map biomes{ {
		{ 0, 0, 204 },
		{ 204, 230, 255 },
		{ 255, 255, 153 },
		{ 170, 255, 128 },
		{ 0, 153, 0 },
		{ 230, 242, 255 },
		} };

	static color_map cold{ {
		{ 230, 247, 255 },
		{ 0, 92, 230 },
		{ 0, 179, 179 },
		{ 0, 51, 153 },
		{ 0, 5, 15 }
		} };

	static color_map warm{ {
		{ 255, 255, 230 },
		{ 255, 204, 0 },
		{ 255, 136, 77 },
		{ 255, 51, 0 },
		{ 128, 0, 0 },
		{ 10, 0, 0 }
		} };

	static color_map quantized{ {
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		}, 6 };

	static color_map pattern{ {
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		{ 255, 255, 255 },
		{ 0, 0, 0 },
		} };

	static const color_map* colormap_presets[] = {
		&jet,
		&classic,
		&grayscale,
		&inv_grayscale,
		&biomes,
		&cold,
		&warm,
		&quantized,
		&pattern
	};

	struct pixel_rgba8 { uint8_t r, g, b, a; };

	inline void make_equalized_histogram(pixel_rgba8* dst, const rs2::video_frame& depth, const color_map& cm)
	{
		SCOPED_PROFILER;

		const size_t max_depth = 0x10000;
		static size_t histogram[max_depth];
		memset(histogram, 0, sizeof(histogram));

		const auto w = depth.get_width(), h = depth.get_height();
		const auto depth_data = (const uint16_t*)depth.get_data();

		for (auto i = 0; i < w*h; ++i) ++histogram[depth_data[i]];
		for (auto i = 2; i < max_depth; ++i) histogram[i] += histogram[i - 1];

		const auto scale = 1.0f / (float)histogram[0xFFFF];

		for (auto i = 0; i < w*h; ++i)
		{
			const auto d = depth_data[i];
			if (d)
			{
				//auto f = histogram[d] / (float)histogram[0xFFFF];
				const auto f = histogram[d] * scale;
				const auto c = cm.get(f);

				dst->r = (uint8_t)c.x;
				dst->g = (uint8_t)c.y;
				dst->b = (uint8_t)c.z;
				dst->a = 255;
			}
			else
			{
				dst->r = 0;
				dst->g = 0;
				dst->b = 0;
				dst->a = 255;
			}
			dst++;
		}
	}

	inline void make_value_cropped_frame(pixel_rgba8* dst, const rs2::video_frame& depth, const color_map& cm, float depth_min, float depth_max, float depth_units)
	{
		SCOPED_PROFILER;

		const auto scale = 1.0f / (depth_max - depth_min);
		const auto w = depth.get_width(), h = depth.get_height();
		const auto depth_data = (const uint16_t*)depth.get_data();

		for (auto i = 0; i < w*h; ++i)
		{
			const auto d = depth_data[i];
			if (d)
			{
				//const auto f = (d * depth_units - depth_min) / (depth_max - depth_min);
				const auto f = (d * depth_units - depth_min) * scale;
				const auto c = cm.get(f);

				dst->r = (uint8_t)c.x;
				dst->g = (uint8_t)c.y;
				dst->b = (uint8_t)c.z;
				dst->a = 255;
			}
			else
			{
				dst->r = 0;
				dst->g = 0;
				dst->b = 0;
				dst->a = 255;
			}
			dst++;
		}
	}

	inline void colorize_depth(pixel_rgba8* dst, const rs2::video_frame& depth, int colormap_id, float depth_min, float depth_max, float depth_units, bool equalize)
	{
		auto& cm = *colormap_presets[colormap_id];

		if (equalize) make_equalized_histogram(dst, depth, cm);
		else make_value_cropped_frame(dst, depth, cm, depth_min, depth_max, depth_units);
	}

};
