// A volumetric rendering of the cornell box

#include "util.frag"

struct Material {
  float3 col;
};

Material RMaterial ( float idx ) {
  if ( idx == 1.0f ) return Material(float3(0.2f, 0.6f, 1.0f));
}

float Line(float3 O, float3 U, float3 V, float R) {
  float3 pa = O - V, ba = U - V;
  float h= clamp(dot(pa, ba)/dot(ba, ba), 0.0, 1.0);
  return length(pa - ba*h) - R;
}

float opQuad ( float3 o, float3 ul, float3 ll, float3 ur, float3 lr) {
  float r = 999.0f;
  opUnion(r, Line(o, ul, ll, 0.005f));
  opUnion(r, Line(o, ll, lr, 0.005f));
  opUnion(r, Line(o, lr, ur, 0.005f));
  opUnion(r, Line(o, ur, ul, 0.005f));
  return r;
}

struct Quad { float3 ul, ll, lr, ur; };
Quad quads[8];

float Quad_Time ( float fi ) {
  return mod(iGlobalTime*2.6f+fi*2.5f, 14.0f)*1.0f - 8.0f;
}

void Initialize_Quad ( ) {
  for ( int i = 0; i != 5; ++ i ) {
    float fi = float(i),
          t  = Quad_Time(fi); 
    float e = Sample_Uniform(fi*fi+fi);
    float3 ul = float3(-1.5f, -1.5f,  t), ll = float3(-1.5f,  1.5f,  t),
           ur = float3( 1.5f, -1.5f,  t), lr = float3( 1.5f,  1.5f,  t);
    opRotate(ul.xy, sin(iGlobalTime)*0.1f);
    opRotate(ul.xz, cos(iGlobalTime*e)*0.2f);
    opRotate(ll.xz, cos(iGlobalTime*e)*0.2f);
    opRotate(lr.xz, cos(iGlobalTime*e)*0.2f);
    opRotate(ur.xz, cos(iGlobalTime*e)*0.2f);

    opRotate(ul.yz, cos(iGlobalTime*e)*0.2f);
    opRotate(ll.yz, cos(iGlobalTime*e)*0.2f);
    opRotate(lr.yz, cos(iGlobalTime*e)*0.2f);
    opRotate(ur.yz, cos(iGlobalTime*e)*0.2f);

    quads[i] = Quad(ul, ll, lr, ur);
  }
}

// Returns <dist, material>
float2 Map ( in vec3 o ) {
  float2 res = float2(999.0f, -1.0f);
  float3 pul, pll, pur, plr;
  for ( int i = 0; i != 5; ++ i ) {
    Quad q = quads[i];
    Union(res, opQuad(o, q.ul, q.ll, q.ur, q.lr), 1.0);
    Quad w = quads[int(mod(float(i+1), 5.0f))];

    float ti = Quad_Time(float(i));
    if ( ti < 2.0f ) {
      Union(res, Line(o, q.ul, w.ul, 0.005f), 1.0);
      Union(res, Line(o, q.ll, w.ll, 0.005f), 1.0);
      Union(res, Line(o, q.ur, w.ur, 0.005f), 1.0);
      Union(res, Line(o, q.lr, w.lr, 0.005f), 1.0);
    }
  }
  return res;
}

float3 BRDF ( float3 O, float3 N, float3 wi, float3 wo, Material mat ) {
  return vec3(clamp(dot(N, wo), 0.0f, 1.0f))*mat.col;
}

void mainImage ( out float4 fragColor, in float2 fragCoord ) {
  float2 uv = -1.0f + 2.0f*fragCoord.xy/iResolution.xy;
  uv.x *= iResolution.x/iResolution.y;
  fragColor = float4(0.0f, 0.0f, 0.0f, 1.0f);

  Initialize_Quad ( );

  float3 col = float3(0.0f);

  Ray eye = Look_At(uv);
  float2 res = March(eye);
  if ( res.x > 0.0f ) {
    float3 O = eye.ori + eye.dir*res.x,
            N = Normal(O),
            wi = eye.dir,
            wo = reflect(wi, N);
    Material mat = RMaterial(res.y);
    col = float3(1.0f);
  }

  fragColor.xyz = col;
}