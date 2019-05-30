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
uniform float iDial6;
uniform float iDial7;
uniform vec2 iResolution;
uniform sampler1D iFFT;
uniform float iNote;
uniform float iPressure;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

void rand(in vec2 x, out float num);
void lfnoise(in vec2 t, out float n);
void stroke(in float d0, in float s, out float d);
void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds);
void rot3(in vec3 p, out mat3 rot);

void dtruchet3(in vec3 x, in float size, out float d)
{
    vec3 y = mod(x, size)-.5*size,
        ind = (x-y)/size;
    
    vec3 r;
    lfnoise(ind.xy-.1*iTime, r.x);
    lfnoise(ind.yz-.1*iTime, r.y);
    lfnoise(ind.zx-.1*iTime, r.z);
    
    mat3 RR;
    rot3(floor(12.*r)*pi/2., RR);
    
    dspline3(RR*y, .5*size*c.zyy, c.yyy, .5*size*c.yzy, d);
    float da;
    dspline3(RR*y, .5*size*c.yyz, c.yyy, .5*size*c.xyy, da);
    d = min(d, da);
    dspline3(RR*y, .5*size*c.yxy, c.yyy, .5*size*c.yyx, da);
    d = min(d, da);
    
    vec3 na, nb;
    lfnoise(62.*x.x*c.xx-1.*iTime, na.x);
    lfnoise(72.*x.y*c.xx-1.3*iTime, na.y);
    lfnoise(-62.*x.x*c.xx+1.*iTime, nb.x);
    lfnoise(-72.*x.y*c.xx+1.3*iTime, nb.y);
    
    stroke(d, .125*size+mix(.005,.007,iScale)*na.x*na.y-mix(.005,01,iScale)*nb.x*nb.y, d);
}

// Scene
float mat;
void scene(in vec3 x, out vec2 d)
{
	x.z -= mix(.1,1.,iDial0)*iTime;
    float n;
    lfnoise(3.*x.z*c.xx-.1*iTime, n);
    n = .5+.5*n;
    
    vec2 cs = vec2(cos(iTime+3.*n), sin(iTime+3.*n));
    mat2 r2 = mat2(cs.x,cs.y,-cs.y,cs.x);
    x.xy = r2 * x.xy;
    
    d = c.xx;
    dtruchet3(x, .1, d.x);
}

// Normal
const float dx = 5.e-4;
void normal(in vec3 x, out vec3 n);

// Texture
void colorize(in vec2 x, out vec3 col)
{    
    float phi = .1*iTime;
    
    vec3 white = vec3(0.89,0.44,0.23),
        gray =vec3(0.25,0.23,0.21);
    float size = .1;
    
    
    vec2 y = mod(x,size)-.5*size;
    y = abs(y)-.001;
    
    float d = min(y.x,y.y);
    col = mix(white, gray, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d*mat));
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
    float d = 0.;//-(o.z-.03)/dir.z;
    int N = 450,
        i;
    
    // Graph
    x = o + d * dir;
    mat3 RR;
            vec3 ra;
            rand(iNote*c.xx+3337., ra.x);
            rand(iNote*c.xx+1337., ra.y);
            rand(iNote*c.xx+2337., ra.z);
            rot3((iNote-48.)*810.*ra,RR);
    // Actual Scene
    {

        // Raymarching
        for(i=0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-4) break;
            d += min(s.x,.01);
        }

        // Illumination
        l = normalize(x-.1*c.yyx);
        if(i<N)
        {
            normal(x,n);
            
            vec3 lo = vec3(0.69,0.23,0.22),
                hi = vec3(1.00,0.33,0.27);
            col  = mix(hi, lo, tanh(d/3.));
            
            col = abs(RR*col);
            col = mix(col, .3*c.xxx, .6);
//             col = mix(col, .3*length(col)/sqrt(3.)*c.xxx, iPressure);
//             col = mix(col, c.yyy, clamp(float(i)/float(N),0.,1.));
        }
    }
    
    // Colorize
    col = .2*col
        + 1.8*col*abs(dot(l,n))
        +1.6*abs(RR*vec3(1.00,0.33,0.27))*(pow(dot(reflect(-l,n),dir),9.));
    col = clamp(col, 0.,1.);
    if(length(col)<.01)col = .1*abs(RR*vec3(0.69,0.23,0.22));
        
    float dd;
    rand(1200.*uv, dd);
    col += dd*.1*c.xxx;
    
    fragColor = clamp(vec4(col,1.0),0.,1.);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
