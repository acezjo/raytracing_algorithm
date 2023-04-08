#include <iostream>
#include <cmath>
#include <cstdlib>

const double pi = 3.1415926535897932385;

class vector {
    public:
        vector() : e{0,0,0} {}
        vector(float e0, float e1, float e2) : e{e0, e1, e2} {}

        float x() const { return e[0]; }
        float y() const { return e[1]; }
        float z() const { return e[2]; }

        vector operator-() const { return vector(-e[0], -e[1], -e[2]); }
        float operator[](int i) const { return e[i]; }
        float& operator[](int i) { return e[i]; }

        vector& operator+=(const vector &v) {
            e[0] += v.e[0];
            e[1] += v.e[1];
            e[2] += v.e[2];
            return *this;
        }

        vector& operator*=(const float t) {
            e[0] *= t;
            e[1] *= t;
            e[2] *= t;
            return *this;
        }

        vector& operator/=(const float t) {
            return *this *= 1/t;
        }

        // vector& operator=(const ppm t) {
        //     e[0] = t.r;
        //     e[1] = t.g;
        //     e[2] = t.b;
        //     return *this;
        // }

        float length() const {
            return sqrt(length_squared());
        }

        float length_squared() const {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }

    public:
        float e[3];
};

inline vector operator+(const vector &u, const vector &v) {
    return vector(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vector operator-(const vector &u, const vector &v) {
    return vector(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vector operator*(const vector &u, const vector &v) {
    return vector(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vector operator*(double t, const vector &v) {
    return vector(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vector operator*(const vector &v, double t) {
    return t * v;
}

inline vector operator/(vector v, double t) {
    return (1/t) * v;
}


inline double dot(const vector &u, const vector &v) {
    return u.e[0] * v.e[0]
         + u.e[1] * v.e[1]
         + u.e[2] * v.e[2];
}

inline vector cross(const vector &u, const vector &v) {
    return vector(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vector normalized(vector v) {
    return v / v.length();
}

inline bool near_zero(const vector& vec) {
    const auto s = 1e-8;
    return (fabs(vec.x()) < s) && (fabs(vec.y()) < s) && (fabs(vec.z()) < s);
}

double random_double(double a, double b) {
    return (double) rand() / RAND_MAX * (b - a) + a;
}

vector random_unit_vector() {
    double a = random_double(0, 2 * pi);
    double z = random_double(-1, 1);
    double r = sqrt(1 - z*z);
    return vector(r*cos(a), r*sin(a), z);
}

vector reflect(const vector& v, const vector& n) {
    return v - 2*dot(v,n)*n;
}

using color = vector;

struct ppm { //Abstracts Pixels with three values into single "objects"
    int r;
    int g;
    int b;

    ppm& operator=(const color t) {
            //maps color vector to gamma corrected 8 bit RGB values for writing to ppm file
            r = 255*sqrt(t.x());
            g = 255*sqrt(t.y());
            b = 255*sqrt(t.z());
            return *this;
        }

    ppm& operator+=(const color t) {
            r = r + 255*t.x();
            g = g + 255*t.y();
            b = b + 255*t.z();
            return *this;
        }

    ppm& operator*=(float t) {
            r = r*t;
            g = g*t;
            b = b*t;
            return *this;
        }
};

struct scattered {
    vector ray_direction;
    color c;
};



struct hit{
    bool is_hit;
    vector intersection_point;
    vector normal_intersection_point;
    color object_color;
};

class Ray {
    public:
        Ray() {}
        Ray(const vector& origin, const vector& direction)
             : orig(origin), dir(direction) 
        {}

        vector origin() const {return orig;}
        vector direction() const {return dir;}

        vector position(float t) const {
            return orig + t*dir;
        }

    public:
        vector orig;
        vector dir;
};



class material {
    public:
    material(){}

    virtual bool emits(const vector& point, color& light){return false;};
    virtual bool scatter(const Ray& ray_in, const hit& intersection, color& attenuation, Ray& scattered) {return false;}};

class diffuse : public material {
    public:
        diffuse(const color& a) : albedo(a){}
        const color albedo;

        virtual bool scatter(const Ray& ray_in, const hit& intersection, color& attenuation, Ray& scattered) {
            vector scatter_direction = intersection.normal_intersection_point + random_unit_vector();
            if (near_zero(scatter_direction)) {
                scatter_direction = intersection.normal_intersection_point;
            }
            scattered = Ray(intersection.intersection_point, scatter_direction);
            attenuation = albedo;
            return true;
        }
};

class specular : public material {
    public:
        specular(const color& r) : reflectance(r){}
        const color reflectance;

        virtual bool scatter(const Ray& r_in, const hit& rec, color& attenuation, Ray& scattered) {
            vector reflected = reflect(normalized(r_in.direction()), rec.normal_intersection_point);
            scattered = Ray(rec.intersection_point, reflected);
            attenuation = reflectance;
            return (dot(scattered.direction(), rec.normal_intersection_point) > 0);
        }
};

class emissive : public material {
public:
    emissive(const color& emitColor, double emitIntensity) :
        m_emitColor(emitColor), m_emitIntensity(emitIntensity) {}

    virtual bool scatter(const Ray& rayIn, const hit& hitRecord, vector& attenuation, Ray& scatteredRay) {
        // Emissive materials don't scatter light, so return false to indicate no further interaction
        return false;
    }

    virtual bool emits(const vector& point, color& light) {
        light = m_emitIntensity * m_emitColor;
        return true;
    }

private:
    color m_emitColor;
    double m_emitIntensity;
};

class camera {
    public:
        camera() {
            aspect_ratio = 16.0 / 9.0;
            viewport_height = 2.0;
            viewport_width = aspect_ratio * viewport_height;
            focal_length = 2.0;

            origin = vector(0, 0, 0);
            horizontal = vector(viewport_width, 0.0, 0.0);
            vertical = vector(0.0, viewport_height, 0.0);
            lower_left_corner = origin - horizontal/2 - vertical/2 - vector(0, 0, focal_length);
        }

        camera(double asp_ratio, double viewp_height, double foc_length) 
             : aspect_ratio(asp_ratio), viewport_height(viewp_height), focal_length(foc_length)
            {
                viewport_width = aspect_ratio * viewport_height;
                origin = vector(0, 0, 0);
                horizontal = vector(viewport_width, 0.0, 0.0);
                vertical = vector(0.0, viewport_height, 0.0);
                lower_left_corner = origin - horizontal/2 - vertical/2 - vector(0, 0, focal_length);
            }

        Ray get_ray(double s, double t) const {
            return Ray(origin, (lower_left_corner + s*horizontal + t*vertical - origin));
        }

    private:
        double aspect_ratio;
        double viewport_height;
        double viewport_width;
        double focal_length;
        vector origin;
        vector lower_left_corner;
        vector horizontal;
        vector vertical;
};

class sphere {
    public:
        sphere(){origin = vector(0,0,-1); r=0.5; c=color(0,0,1);}

        sphere(vector orig, float radius, color color, material* m) 
        : origin(orig), r(radius), c(color), mat(m)
        {}

        hit hit_sphere(Ray ray) {
            hit h;
            vector oc = ray.origin() - origin;
            auto a = dot(ray.direction(), ray.direction());
            auto b = 2.0 * dot(oc, ray.direction());
            auto c = dot(oc, oc) - r*r;
            auto discriminant = b*b - 4*a*c;
            const auto t = (-b - std::sqrt(discriminant)) / (2.0 * a);
            if (discriminant < 0) {
                h.is_hit = false;
            } else {
                if(t < 0) {
                    h.is_hit = false;
                } else {
                    h.intersection_point = ray.origin() + t * ray.direction();
                    h.normal_intersection_point = normalized(h.intersection_point - origin);
                    h.is_hit = true;
                }
            }
            return h;
        }

    public:
        vector origin;
        float r;
        color c;
        material* mat;

};

color ray_color(Ray ray, std::vector<sphere>& spheres, int depth) {
    bool is_hit = false;
    hit closest;
    sphere closest_sphere;
    double closest_distance = std::numeric_limits<double>::max();
    for(auto& sph : spheres) {
        hit h = sph.hit_sphere(ray);
        if(h.is_hit && h.intersection_point.length() < closest_distance) {
            closest = h;
            closest_sphere = sph;
            closest_distance = h.intersection_point.length();
            is_hit = true;
        }
    }
    if (is_hit){
        Ray scattered_ray;
        color attenuation = closest.object_color;
        if (closest_sphere.mat->scatter(ray, closest, attenuation, scattered_ray)) {
            return attenuation * ray_color(scattered_ray, spheres, depth - 1);
        }
        if (closest_sphere.mat->emits(closest.intersection_point, attenuation)) {
            return attenuation;
        }
        return color(0, 0, 0);
    } 
    // code for returning background color:
    vector unit_direction = normalized(ray.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (0.1-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.9, 1.0);

    // code for returning black background
    //return(color(0,0,0));
};