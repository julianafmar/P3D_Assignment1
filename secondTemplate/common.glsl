/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample) //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed); //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    vec3 eyeOffset = cam.eye + cam.u * ls.x + cam.v * ls.y;

    vec3 pixelOnImagePlane = vec3(cam.width * (pixel_sample.x / iResolution.x - 0.5), cam.height * (pixel_sample.y / iResolution.y - 0.5), -cam.planeDist);

    vec3 focalPlaneSample = pixelOnImagePlane * cam.focusDist;

    vec3 rayDirection = normalize(cam.u * (focalPlaneSample.x - ls.x) + cam.v * (focalPlaneSample.y - ls.y) + cam.n * focalPlaneSample.z);

    return createRay(eyeOffset, rayDirection, time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for dialectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDialectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIALECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};


float schlick(float cosine, float refIdx)
{
    float r0 = (1.0 - refIdx) / (1.0 + refIdx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow((1.0 - cosine), 5.0);
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered) {
    if (rec.material.type == MT_DIFFUSE) { // Diffuse reflection
        vec3 target = rec.pos + rec.normal + normalize(randomInUnitSphere(gSeed));
        rScattered = createRay(rec.pos, target - rec.pos);
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;
        return true;
    }

    vec3 reflectNormal;
    if (dot(rIn.d, rec.normal) > 0.0) {
        reflectNormal = rec.normal;
    }
    else {
        reflectNormal = -rec.normal;
    }

    if (rec.material.type == MT_METAL) { // Metal reflection
        vec3 reflected = reflect(normalize(rIn.d), reflectNormal);
        rScattered = createRay(rec.pos, reflected + rec.material.roughness * randomInUnitSphere(gSeed));
        atten = rec.material.specColor;
        return dot(rScattered.d, rec.normal) > 0.0;  
    }
    if (rec.material.type == MT_DIALECTRIC) {
        // Dielectric with refraction and reflection
        atten = vec3(1.0);
        vec3 outwardNormal;        
        float niOverNt;
        vec3 refracted;
        float cosine;
        
        if (dot(rIn.d, rec.normal) > 0.0) { // hit inside
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            cosine = rec.material.refIdx * dot(rIn.d, rec.normal) / length(rIn.d);
            cosine = sqrt(1.0 - rec.material.refIdx * rec.material.refIdx * (1.0 - cosine * cosine));
        } else { // hit from outside
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            cosine = -dot(rIn.d, rec.normal);
        }

        float reflectProb;
        refracted = refract(rIn.d, outwardNormal, niOverNt);
        if (length(refracted) > 0.0) {
            reflectProb = schlick(cosine, rec.material.refIdx);
        } else {
            reflectProb = 1.0;  // Total internal reflection
        }

        if (hash1(gSeed) < reflectProb) {
            vec3 reflected = reflect(rIn.d, rec.normal);
            rScattered = createRay(rec.pos, reflected);
        } else {
            rScattered = createRay(rec.pos, refracted);
        }
        return true;
    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    vec3 v0v1 = t.b - t.a; 
    vec3 v0v2 = t.c - t.a; 
    vec3 normal = normalize(cross(v0v1, v0v2));
    vec3 pvec = cross(r.d, v0v2); 
    float det = dot(v0v1, pvec);

    // Using epsilon to prevent division by zero and culling parallel rays
    if (abs(det) < epsilon) return false;

    float invDet = 1.0 / det; 

    vec3 tvec = r.o - t.a; 
    float u = dot(tvec, pvec) * invDet; 
    if (u < 0.0 || u > 1.0) return false; 

    vec3 qvec = cross(tvec, v0v1); 
    float v = dot(r.d, qvec) * invDet; 
    if (v < 0.0 || u + v > 1.0) return false; 

    float tx = dot(v0v2, qvec) * invDet; 

    if(tx < tmax && tx > tmin)
    {
        rec.t = tx;
        rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}


struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    return mvsphere.center0 + (mvsphere.center1 - mvsphere.center0) * ((time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0));
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec) {
    vec3 oc = s.center - r.o;
    float b = dot(r.d, oc);
    float c = dot(oc, oc) - pow(s.radius, 2.0);

    // Check if ray is outside of the sphere and if sphere is behind the ray, if so return false
    if (c > 0.0 && b <= 0.0) return false; 

    float delta = b * b - c;

    if (delta <= 0.0) return false;

    float t;
    if (c > 0.0) t = b - sqrt(delta);
    else t = b + sqrt(delta);
    
    if (t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, t);
        rec.normal = normalize((rec.pos - s.center) / s.radius);
        return true;
    }
    return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec) {
    float t;
    vec3 center = center(s, r.t);
    vec3 oc = center - r.o;
    float B = dot(r.d, oc);
    float C = dot(oc, oc) - pow(s.radius, 2.0);

    if (C > 0.0 && B <= 0.0) return false;

    float delta = B * B - C;

    if (delta <= 0.0) return false;

    if (C > 0.0) t = B - sqrt(delta);
    else t = B + sqrt(delta);
    
    if (t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, t);
        rec.normal = normalize((rec.pos - center) / s.radius);
        return true;
    }
    return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}
