
#define float4 vec4
#define float3 vec3
#define float2 vec2

#define PI   3.141592653589793f
#define IPI  0.318309886183791f
#define IPI2 0.159154943091895f
#define TAU  6.283185307179586f
#define ITAU 0.159154943091895f

uniform vec3 u_eye3d;
uniform vec3 u_centre3d;
uniform vec3 u_up3d;

struct Ray { float3 ori, dir; };

float2 Map ( float3 o );

float2 March ( in Ray ray ) {
  float dist = 0.0f;
  float2 cur;
  for ( int i = 0; i != 64; ++ i ) {
    cur = Map(ray.ori + ray.dir*dist);
    if ( cur.x <= 0.0001f || dist > 128.0f ) break;
    dist += cur.x;
  }
  if ( dist > 128.0f || dist < 0.0f ) return float2(-1.0f);
  return float2(dist, cur.y);
}

float3 Normal ( float3 p ) {
  float2 e = float2(1.0f, -1.0f)*0.5883f*0.0005f;
  return normalize(
                   e.xyy*Map(p + e.xyy).x +
                   e.yyx*Map(p + e.yyx).x +
                   e.yxy*Map(p + e.yxy).x +
                   e.xxx*Map(p + e.xxx).x);
}

Ray Look_At ( float2 uv ) {
  float3 ori    = float3(-u_eye3d.x,    u_eye3d.y,    u_eye3d.z    ),
         center = float3(-u_centre3d.x, u_centre3d.y, u_centre3d.z ),
         up     = float3(-u_up3d.x,     u_up3d.y,     u_up3d.z     );
  float3 ww = normalize(center - ori),
         uu = normalize(cross(up, ww)),
         vv = normalize(cross(ww, uu));
  return Ray(ori, normalize(uv.x*uu + uv.y*vv + 2.5*ww));
}

// -- maps --
float sdSphere ( float3 o, float r ) { return length(o) - r; }
float sdBox    ( float3 o, float3 b ) {
  float3 d = abs(o) - b;
  return min(max(d.x, max(d.y, d.z)), 0.0f) + length(max(d, 0.0f));
}

void Union ( inout float2 t, float d, in float ID ) {
  if ( t.x > d ) t = float2(d, ID);
}

void opRotate(inout float2 p, in float a ) {
  p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// -- random --
float Sample_Uniform ( float f ) {
  vec2 p = vec2(f, fract(sin(f*23423.234)));
  vec3 p3 = fract(vec3(p.xyx) * 0.1031f);
  p3 += dot(p3, p3.yzx + 19.19f);
  return fract((p3.x + p3.y) * p3.z);
}

float2 Sample_Uniform2 ( float f ) {
  return float2(Sample_Uniform(f),
                Sample_Uniform(Sample_Uniform(f)*23.12310f));
}

// -- sampler --
float3 To_Cartesian ( float cos_theta, float phi ) {
  float sin_theta = sqrt(max(0.0f, 1.0f - cos_theta));
  return float3(cos(phi)*sin_theta, sin(phi)*sin_theta, cos_theta);
}

vec3 Reorient_Hemisphere ( vec3 wo, vec3 N ) {
  vec3 binormal = (abs(N.x) < 1.0 ? vec3(1.0, 0.0, 0.0) :
                   vec3(0.0, 1.0, 0.0));
  binormal = normalize(cross(N, binormal));
  vec3 bitangent = cross(binormal, N);
  return bitangent*wo.x + binormal*wo.y + wo.z*N;
}

vec3 Sample_Cos_Hemisphere ( float r, float3 N, out float pdf ) {
  vec2 u = Sample_Uniform2(r);
  float cos_theta = sqrt(u.y);
  pdf = cos_theta * IPI;
  return Reorient_Hemisphere(normalize(To_Cartesian(cos_theta, TAU*u.x)),
                             N);
}
