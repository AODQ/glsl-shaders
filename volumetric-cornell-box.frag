// A volumetric rendering of the cornell box

#include "util.frag"

// Returns <dist, material>
float2 Map ( in vec3 o ) {
  float2 res = float2(999.0f);
  //-----------wall--------
  // left/right
  Union(res, sdBox(o-vec3( 3.0f,  0.0f,  0.0f), vec3(0.01f, 3.00f, 3.00f)), 1.0f);
  Union(res, sdBox(o-vec3(-3.0f,  0.0f,  0.0f), vec3(0.01f, 3.00f, 3.00f)), 1.0f);
  // up/down
  Union(res, sdBox(o-vec3( 0.0f,  3.0f,  0.0f), vec3(3.00f, 0.01f, 3.00f)), 1.0f);
  Union(res, sdBox(o-vec3( 0.0f, -3.0f,  0.0f), vec3(3.00f, 0.01f, 3.00f)), 1.0f);
  // back
  Union(res, sdBox(o-vec3( 0.0f,  0.0f, -3.0f), vec3(3.00f, 3.00f, 0.01f)), 1.0f);
  //-----------boxes--------
  // -- left box
  vec3 to = o-vec3(1.5f, -1.0f, -0.5f);
  opRotate(to.xz, PI*0.09f);
  Union(res, sdBox(to, vec3(1.0f, 2.0f, 0.6f)), 1.0f);
  // -- right box
  to = o-vec3(-0.8f, -2.0f, 1.5f);
  opRotate(to.xz, -PI*0.19f);
  Union(res, sdBox(to, vec3(1.0f, 1.0f, 0.6f)), 1.0f);
  return res;
}

void mainImage ( out float4 fragColor, in float2 fragCoord ) {
  float2 uv = -1.0f + 2.0f*fragCoord.xy/iResolution.xy;
  uv.x *= iResolution.x/iResolution.y;

  Ray eye = Look_At(uv);
  float2 res = March(eye);

  fragColor = float4(0.0f, 0.0f, 0.0f, 1.0f);
  if ( res.x >= 0.0f ) {
    vec3 Lo = vec3(1.0f);
    vec3 wo = normalize(eye.ori-Lo);
    fragColor.xyz = vec3(clamp(dot(wo, Normal(eye.ori + eye.dir*res.x)), 0.0f, 1.0f));
  }
}
