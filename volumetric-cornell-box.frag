// A volumetric rendering of the cornell box

#include "util.frag"

// Returns <dist, material>
float2 Map ( in vec3 o ) {
  float2 res = float2(999.0f);
  //-----------wall--------
  // left/right
  Union(res, sdBox(o-vec3( 3.0f,  0.0f,  0.0f), vec3(0.01f, 3.00f, 3.00f)), 1.0f);
  Union(res, sdBox(o-vec3(-3.0f,  0.0f,  0.0f), vec3(0.01f, 3.00f, 3.00f)), 2.0f);
  // up/down
  Union(res, sdBox(o-vec3( 0.0f,  3.0f,  0.0f), vec3(3.00f, 0.01f, 3.00f)), 3.0f);
  Union(res, sdBox(o-vec3( 0.0f, -3.0f,  0.0f), vec3(3.00f, 0.01f, 3.00f)), 4.0f);
  // back
  Union(res, sdBox(o-vec3( 0.0f,  0.0f, -3.0f), vec3(3.00f, 3.00f, 0.01f)), 5.0f);
  //-----------boxes--------
  // -- left box
  vec3 to = o-vec3(1.5f, -1.0f, -0.5f);
  opRotate(to.xz, PI*0.09f);
  Union(res, sdBox(to, vec3(1.0f, 2.0f, 0.6f)), 6.0f);
  // -- right box
  to = o-vec3(-0.8f, -2.0f, 1.5f);
  opRotate(to.xz, -PI*0.19f);
  Union(res, sdBox(to, vec3(1.0f, 1.0f, 0.6f)), 7.0f);
  //-----------light--------
  Union(res, sdBox(o-vec3(0.0f, 2.98f, 0.0f), vec3(0.5f, 0.02f, 0.5)), 8.0f);
  //----
  return res;
}

bool Propagate ( inout float3 col, inout float3 final_col, inout Ray eye ) {
  float2 res = March(eye);

  float3 mcol = float3(0.0f);
  if ( res.x < 0.0f ) return false;
  if ( res.y == 1.0f ) mcol = float3(0.3f, 0.8f, 0.3f);
  if ( res.y == 2.0f ) mcol = float3(0.8f, 0.3f, 0.3f);
  if ( res.y == 3.0f ) mcol = float3(0.5f, 0.5f, 0.5f);
  if ( res.y == 4.0f ) mcol = float3(0.5f, 0.5f, 0.5f);
  if ( res.y == 5.0f ) mcol = float3(0.5f, 0.5f, 0.5f);
  if ( res.y == 6.0f ) mcol = float3(0.8f, 0.8f, 0.8f);
  if ( res.y == 7.0f ) mcol = float3(0.6f, 0.6f, 0.6f);

  if ( res.y == 8.0f ) {
    final_col = col;
    return true;
  }

  eye.ori = eye.ori + eye.dir*res.x;
  float3 N = Normal(eye.ori);
  float pdf;
  float3 tdir = Sample_Cos_Hemisphere(dot(eye.ori, eye.dir+iGlobalTime),
                                      Normal(eye.ori), pdf);
  eye.ori += eye.dir*0.001f;

  col = mcol * (abs(dot(eye.dir, N))*IPI)/pdf;
  // -- cast ray temporarily to light source --
  eye.dir = normalize(vec3(0.0f, 2.98f, 0.0f) - eye.ori);
  float3 tcol = mcol * IPI;
  float2 l_res = March(eye);
  if ( l_res.y == 8.0f ) {
    final_col = final_col*0.5f + (col)*0.5f;
  }

  eye.dir = tdir;
  return false;
}

void mainImage ( out float4 fragColor, in float2 fragCoord ) {
  float2 uv = -1.0f + 2.0f*fragCoord.xy/iResolution.xy;
  uv.x *= iResolution.x/iResolution.y;

  Ray eye = Look_At(uv);

  fragColor = float4(0.0f, 0.0f, 0.0f, 1.0f);

  float3 col = float3(1.0f);
  float3 final_col = float3(0.0f);

  for ( int i = 0; i != 5; ++ i ) {
    if ( Propagate(col, final_col, eye) ) {
      final_col /= float(i);
      break;
    }
  }

  fragColor.xyz = (final_col);
}
