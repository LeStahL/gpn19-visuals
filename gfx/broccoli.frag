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
uniform vec2 iResolution;
uniform sampler1D iFFT;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

void rand(in vec2 x, out float num);
void lfnoise(in vec2 t, out float n);

void rand3(in vec3 x, out float num);
void zextrude(in float z, in float d2d, in float h, out float d);
void stroke(in float d0, in float s, out float d);
void smoothmin(in float a, in float b, in float k, out float dst);
void dbox3(in vec3 x, in vec3 b, out float d);

// Random Quadtree
void dcubetree(in vec3 x, in float threshold, in float depth, out float d, out float faco)
{
    d = 1.;
    vec3 y = x, 
        yi;
    float size = .5,
	    fac = 1.;
    faco = 1.;
    for(float i=0.; i<depth; i+=1.)
    {
        vec3 y0 = y;
        y = mod(y, size)-.5*size;
        yi = y0-y;
		float r;
        rand3(yi+fac,r);
        fac *= r*step(r,threshold);
        if(fac != 0.)
        {
            float dd;
            dbox3(y,(.35)*size*c.xxx,dd);
//             dd = mix(dd,length(y)-(.3)*size,step(r,threshold));
            dd = abs(dd)-.01*size;
            smoothmin(d,dd,.01,d);
        } else break;
        
        size *= .5;
    }
    faco += fac*fac;
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
    
    x.z -= .01*iTime;
    
    dcubetree((2.)*x-iTime*c.yyx-.1*iTime, .5,  6.-6.*iScale, d.x, mat);
    float d2, m2 = 0.;
    //dcubetree((3.+.5*iScale)*x+iTime*c.yyx-.1*iTime, .5, 6., d2, m2);
//     mat = 5555.*mat;
	//lfnoise((5.*x.z-5.*iTime)*c.xx, mat);
	//mat = .5+.5*mat;
    //smoothmin(d.x,d2,.2+.2*iScale,d.x);
    
    d2 = length(x.xy)-.1;
    d.x = max(d.x, -d2);
    d = min(d, -length(x.xy)+.3);
    
//     d -= .01;
    
    stroke(d.x, .01, d.x);
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
    float d = (.1)/length(dir.xy);// -(o.z-.12)/dir.z;
    int N = 650,
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
            d += min(s.x,.0005);
        }

        // Illumination
        l = normalize(x+c.yxx);
        if(i<N)
        {
            normal(x,n);
            mat += 14.*sign(n.z);
            col = mix((.5+.5*mat)*c.xxx,(1.+.8*mat)*vec3(0.89,0.44,0.23),.5+.5*sin(x.z));
            col = mix(col,vec3(0.25,0.23,0.21),.5+.5*cos(4.*x.z+mat));
            float phi = atan(x.y, x.x),
                dhex,
                na,
                nal;
            vec2 ind;
            rand(floor(.33*iTime)*c.xx, na);
            rand(floor(.33*iTime)*c.xx+1., nal);
            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));
            
            mat3 RR;
            float ras;
            rand(mat*c.xx,ras);
            rot3(na*1.e3*vec3(1.1,1.5,1.9)+1.*ras+.5*cos(x.z),RR);

            col = mix((.5+.5*mat)*c.xxx,(1.+.8*mat)*abs(RR*vec3(0.89,0.44,0.23)),.5+.5*sin(x.z));
            //rot3(c.xxx+x.z+1200.*na*mat+1.e3*na,RR);
            col = mix(col,abs(RR*vec3(0.25,0.23,0.21)),.5+.5*cos(.5*(x.z)));
            
            col = mix(col, .5*abs(RR*RR*vec3(0.25,0.23,0.21)), clamp(length(x.xy)/.2, 0.,1.));
            col = mix(col, abs(RR*col), step(0.,length(sign(n))));
            
            col = mix(col, 3.*col, step(.9,length(x.xy))*step(1.1,length(x.xy)));
        }
    }
    vec3 c1 = col;
    // Colorize
    col = .8*col
        + .9*col*abs(dot(l,n))
        +5.4*col*abs(pow(dot(reflect(-l,n),dir),3.));
        
    float dd;
    rand(1200.*uv, dd);
    col += dd*.1*c.xxx;
    
    col = mix(col, c1, clamp(float(i)/float(N),0.,1.));
    
    fragColor = clamp(vec4(col,1.0),0.,1.);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
