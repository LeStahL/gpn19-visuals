/* Corfield Imitation 1
 * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

float iScale;

//lengths
float L1 = 2.;
float L2 = 1.;

//masses
float M1 = 3.5;
float M2 = 0.1;

//constants
float G = 9.81;

//runge kutta params
float h = 1.e-1;
float tmax = 3.;

/**eval system of differential equations
params:
tp[0]: theta1
tp[1]: theta2
tp[2]: ptheta1
tp[3]: ptheta2
*/
vec4 f(vec4 tp)
{
    float p0w = (L1*L2*(M1+M2*pow(sin(tp[0]-tp[1]),2.)));
    float C1 = tp[2]*tp[3]*sin(tp[0]-tp[1])/p0w;
    float C2 = (L2*L2*M2*tp[2]*tp[2]+L1*L1*(M1+M2)*tp[3]*tp[3]-L1*L2*M2*tp[2]*tp[3]*cos(tp[0]-tp[1]))*sin(2.*(tp[0]-tp[1]))/(2.*p0w*p0w);
    
    vec4 ret;
    
    ret[0] = (L2*tp[2]-L1*tp[3]*cos(tp[0]-tp[1]))/(L1*p0w);
    ret[1] = (L1*(M1+M2)*tp[3]-L2*M2*tp[2]*cos(tp[0]-tp[1]))/(L2*M2*p0w);
    ret[2] = -(M1+M2)*G*L1*sin(tp[0])-C1+C2;
    ret[3] = -M2*G*L2*sin(tp[1])+C1-C2;
    
    return ret;    
}

vec4 step_rk4(vec4 tp)
{
    vec4 k1 = f(tp);
    vec4 k2 = f(tp + h/2.*k1);
    vec4 k3 = f(tp + h/2.*k2);
    vec4 k4 = f(tp + h*k3);
    return tp + h/6.*(k1+2.*k2+2.*k3+k4);
}

// Hash function
void rand(in vec2 x, out float num)
{
    x += 400.;
    num = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void lfnoise_edge(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    //t = smoothstep(c.yy, c.xx, t);
    t = smoothstep(.4*c.xx,.6*c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}

void rand3(in vec3 x, out float num)
{
    x += 400.;
    num = fract(sin(dot(sign(x)*abs(x) ,vec3(12.9898,78.233,121.112)))*43758.5453);
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// iq's smooth minimum
void smoothmin(in float a, in float b, in float k, out float dst)
{
    float h = max( k-abs(a-b), 0.0 )/k;
    dst = min( a, b ) - h*h*h*k*(1.0/6.0);
}

void dvoronoi(in vec2 x, out float d, out vec2 z)
{
    vec2 y = floor(x);
       float ret = 1.;
    vec2 pf=c.yy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    d = ret;
    z = pf;
}

void dbox3(in vec3 x, in vec3 b, out float d)
{
  vec3 da = abs(x) - b;
  d = length(max(da,0.0))
         + min(max(da.x,max(da.y,da.z)),0.0);
}

void rot3(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

// Scene
float mat;
void scene(in vec3 x, out vec2 d)
{
    d = c.xx;
    
    vec4 state = vec4(x.xy*2.*pi-4.*vec2(pi/4.,pi/4.), 0, 0);
    float time = 0.;
    while(time < tmax) 
    {
        state = step_rk4(state);
        time += h;
    }
    
    d.x = x.z -.02 - .01*log(abs(state.r/state.b));
    //stroke(d.x, .005, d.x);
}

// Normal
const float dx = 5.e-4;
void normal(in vec3 x, out vec3 n)
{
    vec2 s, na;
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

// Texture
void color(in float scale, out vec3 col)
{
    const int N = 5;
    const vec3 colors[N] = vec3[N](
        vec3(1.00,0.67,0.36),
        vec3(0.86,0.45,0.50),
        vec3(0.67,0.43,0.51),
        vec3(0.41,0.36,0.48),
        vec3(0.27,0.36,0.48)
    );
	float index = floor(scale*float(N)), 
        remainder = scale*float(N)-index;
    col = mix(colors[int(index)],colors[int(index)+1], remainder);
}

// Texture
void color2(in float scale, out vec3 col)
{
    const int N = 5;
    const vec3 colors[N] = vec3[N](
        vec3(0.95,0.74,0.56),
        vec3(0.95,0.64,0.50),
        vec3(0.75,0.43,0.42),
        vec3(0.44,0.29,0.36),
        vec3(0.24,0.15,0.23)
    );
	float index = floor(scale*float(N)), 
        remainder = scale*float(N)-index;
    col = mix(colors[int(index)],colors[int(index)+1], remainder);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Set up coordinates
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = c.yyy;
     iScale = mod(iTime,1.);
    uv /= 1.+4.*iScale;
    
    M1 = 2.+1.*sin(iTime);
    M2 = 2.+1.*cos(2.*iTime);
    
    L1 = 2.+.1*sin(iTime);
    L2 = 1.+.1*cos(3.*iTime);
    
    
    /*if(length(uv) > .5)
    {
        fragColor = vec4(col, 0.);
        return;
    }*/
    
   
    
     // Camera setup
    float pp = .3*iTime;
    vec3 o = c.yyx,
        t = c.yyy,
        dir = normalize(t-o),
        r = normalize(c.xyy),
        u = normalize(cross(r,dir)),
        n,
        x,
        l;
    t += uv.x*r + uv.y*u;
    dir = normalize(t-o);
    vec2 s;
    float d = -(o.z-.03)/dir.z;
    int N = 850,
        i;
    
    // Graph
    x = o + d * dir;
    
    // Actual Scene
    {

        // Raymarching
        for(i=0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-4) break;
            d += s.x;//,.005);
        }

        // Illumination
        l = normalize(x+c.yxx);
        if(i<N)
        {
            normal(x,n);
            color((x.z+.01)/.03, col);
            vec3 c1;
            color2((x.z+.01)/.03, c1);
            vec4 state = vec4(x.xy*2.*pi-vec2(pi/2.,pi/2.), 0, 0);
            float time = 0.;
            while(time < tmax) 
            {
                state = step_rk4(state);
                time += h;
            }

            float daa = x.z -.02 - .01*log(abs(state.r/state.b));
            col = mix(col, .7*c1, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, daa));
        }
    }
    
    // Colorize
    col = .2*col
        + 1.3*col*abs(dot(l,n))
        +.4*col*abs(pow(dot(reflect(-l,n),dir),3.));
    
    float dd;
    rand(1200.*uv, dd);
    col += dd*.1*c.xxx;
    
    fragColor = clamp(vec4(col,1.0),0.,1.);
}
