#define PI 3.14159265359
#define AA 3
#define float2 vec2
#define float3 vec3
#define float4 vec4

#define MAX_DIST 1024.0

float max2 ( float2 vec ) {return max(vec.x, vec.y);}
float min2 ( float2 vec ) {return min(vec.x, vec.y);}
float max3 ( float3 vec ) {return max(vec.x, max(vec.y, vec.z));}
float min3 ( float3 vec ) {return min(vec.x, min(vec.y, vec.z));}
float  sqr  ( float  t ) { return t*t; }
float2 sqr2 ( float2 t ) { return t*t; }
float3 sqr3 ( float3 t ) { return t*t; }

float noise ( const in vec3 x ) {
  return sin(dot(x, float3(12.102312, 11.12391, 13.12391))*10000.0);
}

struct Material {
  float metallic, roughness, fresnel, glow, specular;
};

Material RMaterial ( int index ) {
  if ( index == 0 ) return Material (0.1, 0.0, 1.0,  0.1, 1.0); // floor
  if ( index == 1 ) return Material (0.0, 0.0, 1.0,  0.1, 1.0); // wall
  if ( index == 2 ) return Material (1.0, 0.9, 0.09, 0.0, 0.0); // batman
  if ( index == 3 ) return Material (1.0, 0.4, 0.4,  0.0, 0.0); // sphere bump
  if ( index == 4 ) return Material (0.0, 0.1, 1.0,  1.0, 1.0); // specular sphereS
  if ( index == 5 ) return Material (1.0, 0.0, 0.2,  2.0, 0.0); // specular sphereW
  if ( index == 6 ) return Material (1.0, 1.0, 0.0,  1.0, 0.0);
  if ( index == 7 ) return Material (1.0, 1.0, 0.0,  1.0, 0.0);
  if ( index == 8 ) return Material (1.0, 1.0, 0.0,  1.0, 0.0);
  if ( index == 9 ) return Material (1.0, 1.0, 0.0,  1.0, 0.0);
}

float3 RColour ( float3 O, int mat_index ) {
  if ( mat_index == 0 )
    return float3(sin(O.x*PI*0.5)*sin(O.z*PI*0.5) >= 0.0) + float3(0.5);
  if ( mat_index == 1 )
    return float3(sin(O.y*PI*0.5)*sin(O.z*PI*0.5) >= 0.0) + float3(0.5);
  if ( mat_index == 2 )
    return float3(1.0, 1.0, 0.1);
  if ( mat_index == 3 )
    return float3(0.4, 0.9, 0.4);
  if ( mat_index == 4 )
    return float3(0.6);
  if ( mat_index == 5 )
    return float3(0.2);
  return float3(10000.0);
}

const mat2 m2 = mat2(0.6, -0.8, 0.8, 0.6);
const mat3 m3 = mat3(0.00,  0.80,  0.60,
                    -0.80,  0.36, -0.48,
                    -0.60, -0.48,  0.64);
float fbm ( in vec3 p ) {
  float f = 0.0;
  f += 0.5000*noise(p); p = m3*p*2.02;
  f += 0.2500*noise(p); p = m3*p*2.03;
  f += 0.1250*noise(p); p = m3*p*2.01;
  f += 0.0625*noise(p);
  return f/0.9375;
}

float sdPlane ( float3 p, float3 n, float dist ) {
  return dot(p, n) + dist;
}

float sdSphere( float3 p, float s ) {
  return length(p)-s;
}

float sdBox ( float3 p, float3 b ) {
  float3 d = abs(p) - b;
  return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
}


float sdTriPrism ( float3 p, float2 h ) {
  float3 q = abs(p);
  float d1 = q.z-h.y;
  float d2 = max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5;
  return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}


float sdCylinder ( float3 p, float r, float height ) {
  float d = length(p.xz) - r;
  d = max(d, abs(p.y) - height);
  return d;
}

float sdTorus ( float3 p, float hole, float radius ) {
  return length(float2(length(p.xz) - radius, p.y)) - hole;
}

struct MapInfo {
  float dist;
  int mat_index;
  int iters;
};

void opU ( inout MapInfo i, float dist, int mindex ) {
  if ( i.dist > dist ) i = MapInfo(dist, mindex,0);
}

void opUMap ( inout MapInfo i, in MapInfo j ) {
  if ( i.dist > j.dist ) i = j;
}

void opRotate ( inout float2 p, float a ) {
  p = cos(a)*p + sin(a)*float2(p.y, -p.x);
}

MapInfo Streets ( in float3 pos ) {
  MapInfo i = MapInfo(200.0, -1,0);

  // specular sphere 4=s, 5=w
  float A = sdSphere(pos-float3(-1.0,sin(iTime)*1.0+3.0,-4.0), 2.0);
  int mat_index = int(sin(pos.z*PI*1.5)*sin(pos.y*PI*1.5) >= 0.0) + 4;
  opU(i, A, mat_index);

  // bumpy sphere
  pos -= vec3(1.0, 0.0, 5.0);
  opRotate(pos.xy, iTime);
  opRotate(pos.yz, iTime);
  float B = sdSphere(pos, 2.7+sin(pos.x*40.0)*sin(pos.y*40.0)*0.01);
  opU(i, B, 3);
  return i;
}


// :-)
float Batman_Logo ( in float3 p ) {
  MapInfo i = MapInfo(200.0, -1,0);
  float3 g = p;


  { // reflect space
    float3 pnormal = float3(0.0, 0.0, 1.0);
    float t = dot(p, pnormal) + 2.82;
    if ( t < 0.0 ) p = p - (2.0*t)*pnormal;
  }
  opRotate(p.xy, 270.0*PI/180.0);

  float a, b, c, d;

  // ----- wing ----
  a = sdBox(p-vec3( 1.68, 0.0,-1.65),
              vec3(1.7, 0.1, 1.1));
  float h = p.x;
  // -- left half slice
  float3 q = p;
  q.x -= sin(p.y)*0.25;
  q.z += sin(p.x)*0.45;
  b = sdCylinder(q-vec3( 1.2, 0.0,-0.1), 1.2, 0.2);
  // -- bottom half slice
  q = p;
  q.x += sin(p.x)*0.85;
  c = sdCylinder(q-vec3(3.2,0.0,-1.7), 1.2, 0.2);
  float body = max(max(a, -b), -c);
  // wing = c;

  q = p;
  q.x += sin(p.x)*0.3;
  // q.x += cos(p.x)*0.1;
  a = sdCylinder(q-vec3( 0.2, 0.0,-2.5), 0.3, 0.2);
  body = max(body, -a);

  // ------ head ------
  q = p;
  q -= vec3(0.3, 0.0, -2.8);
  opRotate(q.xy, 90.0*PI/180.0);
  opRotate(q.xz, 90.0*PI/180.0);
  opRotate(q.yx, 0.4);
  a = sdTriPrism(q, float2(0.3, 0.1));
  body = min(body, a);
  // body = a;


  // --- bandings ---
  // ring
  q = p;
  a = sdTorus(q-vec3(0.7,0.0,-2.9),0.05, 0.8);
  body = max(max(body, -a), (body-a+0.1)*sqrt(0.5));
  // ridges
  q = p;
  q -= vec3(-0.2, 0.0, -2.0);
  float siz = 1.3;
  q.z = mod(q.z + 0.5*siz, siz) - 0.5*siz;
  a = sdSphere(q, 0.3);
  body = max(max(body, -a), (body-a+0.1)*sqrt(0.5));


  return body;
}

MapInfo Waters ( in float3 pos ) {
  MapInfo i = MapInfo(200.0, -1,0);
  opU(i, pos.y + 1.0, 0);
  opU(i, -pos.y + 39.0, 0);
  opU(i, pos.x + 13.0, 1);
  opU(i, -pos.x + 10.0, 1);
  opU(i, pos.z + 33.0, 0);
  opU(i, -pos.z + 30.0, 0);

  return i;
}

MapInfo Map ( in float3 pos ) {
  MapInfo i = MapInfo(200.0, -1,0);
  opUMap(i, Streets(pos));
  opUMap(i, Waters (pos));

  pos -= float3(-2.0, 3.0, 6.0 + 2.0*cos(iTime));
  float b = Batman_Logo(pos);
  opU(i, b, 2);
  return i;
}

MapInfo March( in float3 ro, in float3 rd ) {
  MapInfo info = MapInfo(0.0, -1,0);
  int i = 0;
  for ( ; i < 256; ++ i ) {
    MapInfo res = Map(ro + rd*info.dist);
    if ( res.dist < 0.01 || info.dist > MAX_DIST ) {
      info.mat_index = res.mat_index;
      break;
    }
    info.dist += res.dist;
  }

  info.iters = i;
  if ( info.dist > MAX_DIST ) info.mat_index = -1;
  return info;
}

mat3 Look_At ( in float3 eye, in float3 center, in float3 up ) {
  float3 cw = normalize(center-eye);
  float3 cu = normalize(cross(cw,up));
  float3 cv = normalize(cross(cu,cw));
  return mat3(cu, cv, cw);
}

float3 Normal ( float3 p ) {
  float2 e = float2(1.0, -1.0)*0.5773*0.01;
  return normalize(
    e.xyy*Map(p + e.xyy).dist +
    e.yyx*Map(p + e.yyx).dist +
    e.yxy*Map(p + e.yxy).dist +
    e.xxx*Map(p + e.xxx).dist);
}

// -----------------------------------------------------------------------------

float3 Illuminate_Norecurse ( float3 O, int ite, float3 Wi, int mat_index );

float3 BRDF_Specular ( float3 origin, int iters, float3 L, float3 V, float3 N, float3 H, Material m,
                       float3 mat_colour) {
  float3 F, G, D;
  // --------- variables

  {// ---------- Fresnel
    // Schlick 1994
    float3 f0 = m.fresnel * (1.0 - m.metallic) + mat_colour*m.metallic;
    F = f0 + (1.0 - f0)*pow(dot(L, H), 5.0);
  }
  {// ---------- Geometry
    // Heits 2014, SmithGGXCorrelated
    float a2 = m.roughness*m.roughness;
    float GGXV = dot(N, L) * sqrt((dot(N, V) - a2*dot(N, V))*dot(N, V) + a2),
          GGXL = dot(N, V) * sqrt((dot(N, L) - a2*dot(N, L))*dot(N, L) + a2);
    G = float3(0.5 / (GGXV + GGXL));
  }
  {// ---------- Distribution
    // Walter et al 2007
    float cosN_H = dot(N, H);
    float a = cosN_H * m.roughness,
          k = m.roughness/(1.0 - (cosN_H*cosN_H));
    D = float3(k*k*(1.0/PI));
  }
    
  float3 val = float3(F*G*D);
  if ( m.specular > 0.0 ) {
    float3 R = normalize(reflect(-V, N));
    origin += R*1.0;
    MapInfo info = March(origin, R);
    val = Illuminate_Norecurse(origin + R*info.dist, iters, R, info.mat_index);
  }

  // ---------- add it all up

  return val;
}

float3 BRDF_Diffuse ( float3 L, float3 V, float3 N, float3 H, Material m,
                      float3 mat_colour) {
  // Burley 2012, Physically-Based shading at disney
  float f90 = 0.5 + 2.0*m.roughness * sqr(dot(L, H));
  float light_scatter = 1.0 + (f90 - 1.0) * pow(1.0 - dot(L, N), 5.0),
        view_scatter  = 1.0 + (f90 - 1.0) * pow(1.0 - dot(V, N), 5.0);
  float3 D = float3(light_scatter * view_scatter * (1.0/PI));
  return D * (1.0 - m.metallic) * mat_colour;
}

float Shadow ( in float3 ro, in vec3 rd, vec3 lo ) {
  float max_t = length(ro - lo);
  float res = 1.0;
  float t = 0.02;
  for ( int i = 0; i != 24; ++ i ) {
    float dist = Map(ro + rd*t).dist;
    if ( dist < 0.001 )
      return 0.0;
    t += dist;
    res = min(res, 512.0 * dist/t);
    if ( t >= 20.5 ) break;
  }
  return clamp(res, 0.0, 1.0);
}

float3 Illuminate_Norecurse ( float3 O, int iters, float3 Wi, int mat_index ) {
  float3 colour = float3(0.0);
  for ( int i = 0; i != 2; ++ i ) {
    // ------- lights
    float3 lo, emission;
    if ( i == 0 ) {
      float tim = iTime*0.3;
      lo = float3(cos(iTime)*6.0 - 2.0, 7.0, sin(iTime)*-5.0);
      emission = float3(0.7);
    }
    if ( i == 1 ) {
      lo = float3(cos(iTime) - 2.0, 7.0, 12.0 - cos(iTime)*3.0);
      emission = float3(0.4, 0.03, 0.03);
    }
    // ------- vars
    float3 V = normalize(-Wi),
           N = Normal(O),
           L = normalize(lo - O),
           H = normalize(V + L),
           R = normalize(reflect(Wi, N));
    Material material = RMaterial(mat_index);
    float3 mat_colour = RColour(O, mat_index);
    // --- specular BRDF ---
    // float3 brdf_sp = BRDF_Specular(O, iters, L, V, N, H, material, mat_colour);
    // --- diffuse  BRDF ---
    float3 brdf_di = BRDF_Diffuse(L, V, N, H, material, mat_colour);

    colour += float(iters) * 0.005 * material.glow;
    colour += (brdf_di) * dot(N, L) * emission;
  }
  return colour/2.0;
}

float3 Illuminate ( float3 O, int iters, float3 Wi, int mat_index ) {
  float3 colour = float3(0.0);
  for ( int i = 0; i != 2; ++ i ) {
    // ------- lights
    float3 lo, emission;
    if ( i == 0 ) {
      float tim = iTime*0.3;
      lo = float3(cos(iTime)*5.0, 7.0, sin(iTime)*-5.0);
      emission = float3(0.7);
    }
    if ( i == 1 ) {
      lo = float3(cos(iTime), 7.0, 12.0 - cos(iTime)*3.0);
      emission = float3(0.4, 0.03, 0.03);
    }
    // ------- vars
    float3 V = normalize(-Wi),
           N = Normal(O),
           L = normalize(lo - O),
           H = normalize(V + L),
           R = normalize(reflect(Wi, N));
    Material material = RMaterial(mat_index);
    float3 mat_colour = RColour(O, mat_index);
    // --- specular BRDF ---
    float3 brdf_sp = BRDF_Specular(O, iters, L, V, N, H, material, mat_colour);
    // --- diffuse  BRDF ---
    float3 brdf_di = BRDF_Diffuse(L, V, N, H, material, mat_colour);

    colour += float(iters) * 0.005 * material.glow;

    colour += (brdf_di + brdf_sp) * dot(N, L) * emission;
  }
  return colour/2.0;
}
// ----------------------------------------------------------------------------
void mainImage( out float4 frag_colour, in float2 frag_coord ) {
  float3 colour = float3(0.0);
  for ( int i = 0; i != AA; ++ i )
  for ( int j = 0; j != AA; ++ j ) {
    float2 rr = float2(float(i), float(j))/float(AA);

    float2 pixel = (-iResolution.xy + 2.0*(frag_coord+rr))/iResolution.y;
    float3 ro = float3(8.0, 6.0 + cos(iTime)*4.0,
                      1.0 + (sin(iTime*0.5))*8.0),
          rd = Look_At(ro, float3(0.0), float3(0.0, 1.0, 0.0))
                *normalize(float3(pixel, 2.0));
    MapInfo info = March(ro, rd);
    if ( info.mat_index < 0 ) {
      frag_colour = vec4(0.0);
    }
    else {
      colour += Illuminate(ro + rd*info.dist, info.iters, rd, info.mat_index);
    }
  }
  frag_colour = vec4(colour/float(AA*AA), 1.0);
}

