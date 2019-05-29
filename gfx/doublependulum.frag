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
  
#version 130

uniform float iTime;
uniform float iFFTWidth;
uniform float iScale;
uniform float iHighScale;
uniform float iNBeats;
uniform float iDial0;
uniform float iDial7;
uniform vec2 iResolution;
uniform sampler1D iFFT;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

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

void rand(in vec2 x, out float num);
void lfnoise(in vec2 t, out float n);
void stroke(in float d0, in float s, out float d);
void rot3(in vec3 p, out mat3 rot);

// Scene
float mat;
void scene(in vec3 x, out vec2 d)
{
    d = c.xx;
    
    x.y -= .1*iTime;
    
    vec4 state = vec4(x.xy*2.*pi-vec2(pi,pi), 0, 0);
    float time = 0.;
    while(time < tmax) 
    {
        state = step_rk4(state);
        time += h;
    }
    
    d.x = x.z -.02 + .01*abs(log(abs(state.r/state.b)));
    //stroke(d.x, .005, d.x);
}

// Normal
const float dx = 5.e-4;
void normal(in vec3 x, out vec3 n);

// Texture
void color(in float scale, out vec3 col)
{
    scale = fract(scale);
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
    col = fract(col);
}

// Texture
void color2(in float scale, out vec3 col)
{
    scale = fract(scale);
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
    col = fract(col);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Set up coordinates
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = c.yyy;

    if(length(uv) > .5)
    {
        fragColor = vec4(col, 0.);
        return;
    }
    vec3 na;
    lfnoise(iTime*c.xx, na.x);
    lfnoise(2.*iTime*c.xx+2337., na.y);
    lfnoise(3.*iTime*c.xx+3337., na.z);
    
    M1 = 2.+mix(0.,1.,iScale)*na.x;
    M2 = 2.+mix(0.,1.,iScale)*na.y;
    
    L1 = 2.+.1*mix(0.,1.,iScale)*na.x;
    L2 = 1.+.1*mix(0.,1.,iScale)*na.z;
    
    
    uv /= mix(.5,8.,iDial0);
    
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
    float d = -(o.z)/dir.z;
    int N = 450,
        i;
    
    // Graph
    x = o + d * dir;
    
    // Actual Scene
    {

        // Raymarching
//         for(i=0; i<N; ++i)
//         {
//             x = o + d * dir;
//             scene(x,s);
//             if(s.x < 1.e-4) break;
//             d += s.x;//,.005);
//         }

        // Illumination
        l = normalize(x+c.yxx);
        if(i<N)
        {
            vec4 state = vec4((x.xy-.1*iTime*c.yx)*2.*pi-vec2(pi,pi), 0, 0);
            float time = 0.;
            while(time < tmax) 
            {
                state = step_rk4(state);
                time += h;
            }
            float da = -.02 - .01*log(abs(state.r/state.b));
            d += da;
            x = o + d * dir;
            normal(x,n);
            color(abs(x.z)*4., col);
            vec3 c1 = c.yyy;
            color2(abs(x.z)*4., c1);
            
            col = mix(col, c1, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, abs(da)-.1));
            
                        float na, nal;
			rand(floor(.33*iTime)*c.xx, na);
            rand(floor(.33*iTime)*c.xx+1., nal);
            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));
             mat3 RR;
            rot3(na*1.e3*vec3(1.1,1.5,1.9)+13.*length(col),RR);
            col = abs(RR*col);
            
            col = mix(col, length(col)/sqrt(3)*c.xxx, .7);
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

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
