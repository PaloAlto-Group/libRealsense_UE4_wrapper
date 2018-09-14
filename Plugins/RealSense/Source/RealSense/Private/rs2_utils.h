// License: Apache 2.0. See LICENSE file in root directory.
// Copyright(c) 2015 Intel Corporation. All Rights Reserved.

#pragma once

#if !defined(M_PI)
#define M_PI 3.14159265358979323846264338327950288
#endif

namespace rs2
{
    inline float clamp(float x, float min, float max)
    {
        return std::max(std::min(max, x), min);
    }

    inline float smoothstep(float x, float min, float max)
    {
        x = clamp((x - min) / (max - min), 0.0, 1.0);
        return x*x*(3 - 2 * x);
    }

    inline float lerp(float a, float b, float t)
    {
        return b * t + a * (1 - t);
    }

    struct plane
    {
        float a;
        float b;
        float c;
        float d;
    };
    inline bool operator==(const plane& lhs, const plane& rhs) { return lhs.a == rhs.a && lhs.b == rhs.b && lhs.c == rhs.c && lhs.d == rhs.d; }

    struct float3
    {
        float x, y, z;

        float length() const { return sqrt(x*x + y*y + z*z); }

        float3 normalize() const
        {
            return (length() > 0)? float3{ x / length(), y / length(), z / length() }:*this;
        }
    };

    inline float3 cross(const float3& a, const float3& b)
    {
        return { a.y * b.z - b.y * a.z, a.x * b.z - b.x * a.z, a.x * b.y - a.y * b.x };
    }

    inline float evaluate_plane(const plane& plane, const float3& point)
    {
        return plane.a * point.x + plane.b * point.y + plane.c * point.z + plane.d;
    }

    inline float3 operator*(const float3& a, float t)
    {
        return { a.x * t, a.y * t, a.z * t };
    }

    inline float3 operator/(const float3& a, float t)
    {
        return { a.x / t, a.y / t, a.z / t };
    }

    inline float3 operator+(const float3& a, const float3& b)
    {
        return { a.x + b.x, a.y + b.y, a.z + b.z };
    }

    inline float3 operator-(const float3& a, const float3& b)
    {
        return { a.x - b.x, a.y - b.y, a.z - b.z };
    }

    inline float3 lerp(const float3& a, const float3& b, float t)
    {
        return b * t + a * (1 - t);
    }

    struct float2
    {
        float x, y;

        float length() const { return sqrt(x*x + y*y); }

        float2 normalize() const
        {
            return { x / length(), y / length() };
        }
    };

    inline float3 lerp(const std::array<float3, 4>& rect, const float2& p)
    {
        auto v1 = lerp(rect[0], rect[1], p.x);
        auto v2 = lerp(rect[3], rect[2], p.x);
        return lerp(v1, v2, p.y);
    }

    using plane_3d = std::array<float3, 4>;

    inline std::vector<plane_3d> subdivide(const plane_3d& rect, int parts = 4)
    {
        std::vector<plane_3d> res;
        res.reserve(parts*parts);
        for (float i = 0.f; i < parts; i++)
        {
            for (float j = 0.f; j < parts; j++)
            {
                plane_3d r;
                r[0] = lerp(rect, { i / parts, j / parts });
                r[1] = lerp(rect, { i / parts, (j + 1) / parts });
                r[2] = lerp(rect, { (i + 1) / parts, (j + 1) / parts });
                r[3] = lerp(rect, { (i + 1) / parts, j / parts });
                res.push_back(r);
            }
        }
        return res;
    }

    inline float operator*(const float3& a, const float3& b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    inline bool is_valid(const plane_3d& p)
    {
        std::vector<float> angles;
        angles.reserve(4);
        for (int i = 0; i < p.size(); i++)
        {
            auto p1 = p[i];
            auto p2 = p[(i+1) % p.size()];
            if ((p2 - p1).length() < 1e-3) return false;

            p1 = p1.normalize();
            p2 = p2.normalize();

            angles.push_back(acos((p1 * p2) / sqrt(p1.length() * p2.length())));
        }
        return std::all_of(angles.begin(), angles.end(), [](float f) { return f > 0; }) ||
               std::all_of(angles.begin(), angles.end(), [](float f) { return f < 0; });
    }

    inline float2 operator-(float2 a, float2 b)
    {
        return { a.x - b.x, a.y - b.y };
    }

    inline float2 operator*(float a, float2 b)
    {
        return { a * b.x, a * b.y };
    }

    struct matrix4
    {
        float mat[4][4];

        matrix4()
        {
            std::memset(mat, 0, sizeof(mat));
        }

        matrix4(float vals[4][4])
        {
            std::memcpy(mat,vals,sizeof(mat));
        }

        //init rotation matrix from quaternion
        matrix4(rs2_quaternion q)
        {
            mat[0][0] = 1 - 2*q.y*q.y - 2*q.z*q.z; mat[0][1] = 2*q.x*q.y - 2*q.z*q.w;     mat[0][2] = 2*q.x*q.z + 2*q.y*q.w;     mat[0][3] = 0.0f;
            mat[1][0] = 2*q.x*q.y + 2*q.z*q.w;     mat[1][1] = 1 - 2*q.x*q.x - 2*q.z*q.z; mat[1][2] = 2*q.y*q.z - 2*q.x*q.w;     mat[1][3] = 0.0f;
            mat[2][0] = 2*q.x*q.z - 2*q.y*q.w;     mat[2][1] = 2*q.y*q.z + 2*q.x*q.w;     mat[2][2] = 1 - 2*q.x*q.x - 2*q.y*q.y; mat[2][3] = 0.0f;
            mat[3][0] = 0.0f;                      mat[3][1] = 0.0f;                      mat[3][2] = 0.0f;                      mat[3][3] = 1.0f;
        }

        //init translation matrix from vector
        matrix4(rs2_vector t)
        {
            mat[0][0] = 1.0f; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = t.x;
            mat[1][0] = 0.0f; mat[1][1] = 1.0f; mat[1][2] = 0.0f; mat[1][3] = t.y;
            mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = 1.0f; mat[2][3] = t.z;
            mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;
        }

        rs2_quaternion normalize(rs2_quaternion a)
        {
            float norm = sqrtf(a.x*a.x + a.y*a.y + a.z*a.z + a.w*a.w);
            rs2_quaternion res = a;
            res.x /= norm;
            res.y /= norm;
            res.z /= norm;
            res.w /= norm;
            return res;
        }

        rs2_quaternion to_quaternion()
        {
            float tr[4];
            rs2_quaternion res;
            tr[0] = (mat[0][0] + mat[1][1] + mat[2][2]);
            tr[1] = (mat[0][0] - mat[1][1] - mat[2][2]);
            tr[2] = (-mat[0][0] + mat[1][1] - mat[2][2]);
            tr[3] = (-mat[0][0] - mat[1][1] + mat[2][2]);
            if (tr[0] >= tr[1] && tr[0] >= tr[2] && tr[0] >= tr[3])
            {
                float s = 2 * sqrt(tr[0] + 1);
                res.w = s / 4;
                res.x = (mat[2][1] - mat[1][2]) / s;
                res.y = (mat[0][2] - mat[2][0]) / s;
                res.z = (mat[1][0] - mat[0][1]) / s;
            }
            else if (tr[1] >= tr[2] && tr[1] >= tr[3]) {
                float s = 2 * sqrt(tr[1] + 1);
                res.w = (mat[2][1] - mat[1][2]) / s;
                res.x = s / 4;
                res.y = (mat[1][0] + mat[0][1]) / s;
                res.z = (mat[2][0] + mat[0][2]) / s;
            }
            else if (tr[2] >= tr[3]) {
                float s = 2 * sqrt(tr[2] + 1);
                res.w = (mat[0][2] - mat[2][0]) / s;
                res.x = (mat[1][0] + mat[0][1]) / s;
                res.y = s / 4;
                res.z = (mat[1][2] + mat[2][1]) / s;
            }
            else {
                float s = 2 * sqrt(tr[3] + 1);
                res.w = (mat[1][0] - mat[0][1]) / s;
                res.x = (mat[0][2] + mat[2][0]) / s;
                res.y = (mat[1][2] + mat[2][1]) / s;
                res.z = s / 4;
            }
            return normalize(res);
        }

        void to_column_major(float column_major[16])
        {
            column_major[0] = mat[0][0];
            column_major[1] = mat[1][0];
            column_major[2] = mat[2][0];
            column_major[3] = mat[3][0];
            column_major[4] = mat[0][1];
            column_major[5] = mat[1][1];
            column_major[6] = mat[2][1];
            column_major[7] = mat[3][1];
            column_major[8] = mat[0][2];
            column_major[9] = mat[1][2];
            column_major[10] = mat[2][2];
            column_major[11] = mat[3][2];
            column_major[12] = mat[0][3];
            column_major[13] = mat[1][3];
            column_major[14] = mat[2][3];
            column_major[15] = mat[3][3];
        }
    };

    inline matrix4 operator*(matrix4 a, matrix4 b)
    {
        matrix4 res;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                float sum = 0.0f;
                for (int k = 0; k < 4; k++)
                {
                    sum += a.mat[i][k] * b.mat[k][j];
                }
                res.mat[i][j] = sum;
            }
        }
        return res;
    }

    inline matrix4 tm2_pose_to_world_transformation(const rs2_pose& pose)
    {
        matrix4 rotation(pose.rotation);
        matrix4 translation(pose.translation);
        matrix4 G_tm2_body_to_tm2_world = translation * rotation;
        float rotate_180_y[4][4] = { { -1, 0, 0, 0 },
                                     { 0, 1, 0, 0 },
                                     { 0, 0,-1, 0 },
                                     { 0, 0, 0, 1 } };
        matrix4 G_vr_body_to_tm2_body(rotate_180_y);
        matrix4 G_vr_body_to_tm2_world = G_tm2_body_to_tm2_world * G_vr_body_to_tm2_body;

        float rotate_90_x[4][4] = { { 1, 0, 0, 0 },
                                    { 0, 0,-1, 0 },
                                    { 0, 1, 0, 0 },
                                    { 0, 0, 0, 1 } };
        matrix4 G_tm2_world_to_vr_world(rotate_90_x);
        matrix4 G_vr_body_to_vr_world = G_tm2_world_to_vr_world * G_vr_body_to_tm2_world;

        return G_vr_body_to_vr_world;
    }

    inline rs2_pose correct_tm2_pose(const rs2_pose& pose)
    {
        matrix4 G_vr_body_to_vr_world = tm2_pose_to_world_transformation(pose);
        rs2_pose res = pose;
        res.translation.x = G_vr_body_to_vr_world.mat[0][3];
        res.translation.y = G_vr_body_to_vr_world.mat[1][3];
        res.translation.z = G_vr_body_to_vr_world.mat[2][3];
        res.rotation = G_vr_body_to_vr_world.to_quaternion();
        return res;
    }

    // return the distance between p and the line created by p1 and p2
    inline float point_to_line_dist(float2 p1, float2 p2, float2 p)
    {
        float d = abs((p2.x - p1.x)*(p1.y - p.y) - (p1.x - p.x)*(p2.y - p1.y)) / sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
        return d;
    }

    inline std::vector<float2> simplify_line(const std::vector<float2>& points)
    {
        std::vector<float2> res;
        float max_distance = 0.0f;
        int max_distance_index = 0;
        float distance_limit = 0.01f; //1 centimeter
        // Find the point with the maximum distance from the 2 end points of the vector
        for (int i = 1; i < points.size() - 1; i++)
        {
            float d = point_to_line_dist(points[0], points.back(), points[i]);
            if (d > max_distance)
            {
                max_distance = d;
                max_distance_index = i;
            }
        }
        // If max distance is greater than the limit, recursively simplify
        if (max_distance > distance_limit)
        {
            // Recursive call
            std::vector<float2> first_half(points.begin(), points.begin() + max_distance_index);
            std::vector<float2> second_half(points.begin() + max_distance_index, points.end());
            res = simplify_line(first_half);
            std::vector<float2> res_second_half = simplify_line(second_half);
            //check if the connection points of the 2 halves are too close
            float2 p1 = res.back();
            float2 p2 = res_second_half[0];
            if (sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2)) < 0.01)
            {
                res.insert(res.end(), res_second_half.begin() + 1, res_second_half.end());
            }
            else
            {
                res.insert(res.end(), res_second_half.begin(), res_second_half.end());
            }
        }
        else
        {
            res.push_back(points[0]);
            res.push_back(points.back());
        }

        return res;
    }

    inline bool point_in_polygon_2D(const std::vector<float2>& polygon, float2 point)
    {
        bool inside = false;
        int i = 0, j = 0;
        for (i = 0, j = static_cast<int>(polygon.size()) - 1; i < static_cast<int>(polygon.size()); j = i++)
        {
            if (((polygon[i].y > point.y) != (polygon[j].y > point.y)) &&
                (point.x < (polygon[j].x - polygon[i].x) * (point.y - polygon[i].y) / (polygon[j].y - polygon[i].y) + polygon[i].x))
            {
                inside = !inside;
            }
        }
        return inside;
    }

    template<typename T>
    T normalizeT(const T& in_val, const T& min, const T& max)
    {
        return ((in_val - min)/(max - min));
    }

    template<typename T>
    T unnormalizeT(const T& in_val, const T& min, const T& max)
    {
        if (min == max) return min;
        return ((in_val * (max - min)) + min);
    }

    struct rect
    {
        float x, y;
        float w, h;

        void operator=(const rect& other)
        {
            x = other.x;
            y = other.y;
            w = other.w;
            h = other.h;
        }

        operator bool() const
        {
            return w*w > 0 && h*h > 0;
        }

        bool operator==(const rect& other) const
        {
            return x == other.x && y == other.y && w == other.w && h == other.h;
        }

        bool operator!=(const rect& other) const
        {
            return !(*this == other);
        }

        rect normalize(const rect& normalize_to) const
        {
            return rect{normalizeT(x, normalize_to.x, normalize_to.x + normalize_to.w),
                        normalizeT(y, normalize_to.y, normalize_to.y + normalize_to.h),
                        normalizeT(w, 0.f, normalize_to.w),
                        normalizeT(h, 0.f, normalize_to.h)};
        }

        rect unnormalize(const rect& unnormalize_to) const
        {
            return rect{unnormalizeT(x, unnormalize_to.x, unnormalize_to.x + unnormalize_to.w),
                        unnormalizeT(y, unnormalize_to.y, unnormalize_to.y + unnormalize_to.h),
                        unnormalizeT(w, 0.f, unnormalize_to.w),
                        unnormalizeT(h, 0.f, unnormalize_to.h)};
        }

        // Calculate the intersection between two rects
        // If the intersection is empty, a rect with width and height zero will be returned
        rect intersection(const rect& other) const
        {
            auto x1 = std::max(x, other.x);
            auto y1 = std::max(y, other.y);
            auto x2 = std::min(x + w, other.x + other.w);
            auto y2 = std::min(y + h, other.y + other.h);

            return{
                x1, y1,
                std::max(x2 - x1, 0.f),
                std::max(y2 - y1, 0.f)
            };
        }

        // Calculate the area of the rect
        float area() const
        {
            return w * h;
        }

        rect cut_by(const rect& r) const
        {
            auto x1 = x;
            auto y1 = y;
            auto x2 = x + w;
            auto y2 = y + h;

            x1 = std::max(x1, r.x);
            x1 = std::min(x1, r.x + r.w);
            y1 = std::max(y1, r.y);
            y1 = std::min(y1, r.y + r.h);

            x2 = std::max(x2, r.x);
            x2 = std::min(x2, r.x + r.w);
            y2 = std::max(y2, r.y);
            y2 = std::min(y2, r.y + r.h);

            return { x1, y1, x2 - x1, y2 - y1 };
        }

        bool contains(const float2& p) const
        {
            return (p.x >= x) && (p.x < x + w) && (p.y >= y) && (p.y < y + h);
        }

        rect pan(const float2& p) const
        {
            return { x - p.x, y - p.y, w, h };
        }

        rect center() const
        {
            return{ x + w / 2.f, y + h / 2.f, 0, 0 };
        }

        rect lerp(float t, const rect& other) const
        {
            return{
                rs2::lerp(x, other.x, t), rs2::lerp(y, other.y, t),
                rs2::lerp(w, other.w, t), rs2::lerp(h, other.h, t),
            };
        }

        rect adjust_ratio(float2 size) const
        {
            auto H = static_cast<float>(h), W = static_cast<float>(h) * size.x / size.y;
            if (W > w)
            {
                auto scale = w / W;
                W *= scale;
                H *= scale;
            }

            return{ float(floor(x + floor(w - W) / 2)),
                    float(floor(y + floor(h - H) / 2)),
                    W, H };
        }

        rect scale(float factor) const
        {
            return { x, y, w * factor, h * factor };
        }

        rect grow(int pixels) const
        {
            return { x - pixels, y - pixels, w + pixels*2, h + pixels*2 };
        }

        rect grow(int dx, int dy) const
        {
            return { x - dx, y - dy, w + dx*2, h + dy*2 };
        }

        rect shrink_by(float2 pixels) const
        {
            return { x + pixels.x, y + pixels.y, w - pixels.x * 2, h - pixels.y * 2 };
        }

        rect center_at(const float2& new_center) const
        {
            auto c = center();
            auto diff_x = new_center.x - c.x;
            auto diff_y = new_center.y - c.y;

            return { x + diff_x, y + diff_y, w, h };
        }

        rect fit(rect r) const
        {
            float new_w = w;
            float new_h = h;

            if (w < r.w)
                new_w = r.w;

            if (h < r.h)
                new_h = r.h;

            auto res = rect{x, y, new_w, new_h};
            return res.adjust_ratio({w,h});
        }

        rect zoom(float zoom_factor) const
        {
            auto c = center();
            return scale(zoom_factor).center_at({c.x,c.y});
        }

        rect enclose_in(rect in_rect) const
        {
            rect out_rect{x, y, w, h};
            if (w > in_rect.w || h > in_rect.h)
            {
                return in_rect;
            }

            if (x < in_rect.x)
            {
                out_rect.x = in_rect.x;
            }

            if (y < in_rect.y)
            {
                out_rect.y = in_rect.y;
            }


            if (x + w > in_rect.x + in_rect.w)
            {
                out_rect.x = in_rect.x + in_rect.w - w;
            }

            if (y + h > in_rect.y + in_rect.h)
            {
                out_rect.y = in_rect.y + in_rect.h - h;
            }

            return out_rect;
        }

        bool intersects(const rect& other) const
        {
            return other.contains({ x, y }) || other.contains({ x + w, y }) ||
                other.contains({ x, y + h }) || other.contains({ x + w, y + h }) ||
                contains({ other.x, other.y });
        }
    };

    inline bool is_integer(float f)
    {
        return (fabs(fmod(f, 1)) < std::numeric_limits<float>::min());
    }

    // Comparing parameter against a range of values of the same type
    // https://stackoverflow.com/questions/15181579/c-most-efficient-way-to-compare-a-variable-to-multiple-values
    template <typename T>
    bool val_in_range(const T& val, const std::initializer_list<T>& list)
    {
        for (const auto& i : list) {
            if (val == i) {
                return true;
            }
        }
        return false;
    }

    // RS4xx with RealTec RGB sensor may additionally require sensor orientation control to make runtime adjustments
    inline void rotate_rgb_image(device& dev,uint32_t res_width)
    {
        static bool flip = true;
        uint8_t hor_flip_val{}, ver_flip_val{};

        if (flip)
        {
            hor_flip_val = ((res_width < 1280) ? (uint8_t)0x84 : (uint8_t)0x20);
            ver_flip_val = ((res_width < 1280) ? (uint8_t)0x47 : (uint8_t)0x46);
        }
        else
        {
            hor_flip_val = ((res_width < 1280) ? (uint8_t)0x82 : (uint8_t)0x86);
            ver_flip_val = ((res_width < 1280) ? (uint8_t)0x41 : (uint8_t)0x40);
        }

        std::vector<uint8_t> hor_flip{ 0x14, 0, 0xab, 0xcd, 0x29, 0, 0, 0, 0x20, 0x38, 0x0, 0x0,
            hor_flip_val, 0,0,0,0,0,0,0,0,0,0,0 };
        std::vector<uint8_t> ver_flip{ 0x14, 0, 0xab, 0xcd, 0x29, 0, 0, 0, 0x21, 0x38, 0x0, 0x0,
            ver_flip_val, 0,0,0,0,0,0,0,0,0,0,0 };

        dev.as<debug_protocol>().send_and_receive_raw_data(hor_flip);
        dev.as<debug_protocol>().send_and_receive_raw_data(ver_flip);

        flip = !flip;
    }
}
