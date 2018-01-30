
#define float4 vec4
#define float3 vec3
#define float2 vec2

#define PI  (3.141592654f)
#define TAU (6.283185307f)

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