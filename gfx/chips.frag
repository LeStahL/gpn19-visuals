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

void rand(in vec2 x, out float num);
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

void zextrude(in float z, in float d2d, in float h, out float d);
void stroke(in float d0, in float s, out float d);
void smoothmin(in float a, in float b, in float k, out float dst);
void dvoronoi(in vec2 x, out float d, out vec2 z);
void rot3(in vec3 p, out mat3 rot);

float mat;
void scene(in vec3 x, out vec2 d)
{
	x.y += .1*mix(1.,5.,iDial0)*iTime;
    
    d = c.xx;
    d.x = x.z;
    
	vec3 n;
    lfnoise_edge(34.*x.xy, n.x);
    lfnoise_edge(14.*x.xy+1337., n.y);
    lfnoise_edge(54.*x.xy+2337., n.z);
    
    float da;
    stroke(2.*n.x-3.*n.y-n.z, .02, da);
    
    zextrude(x.z, -da, .02, da);
    stroke(da, mix(.01,.05,iScale), da);
    stroke(da,mix(.1,.5,iScale),da);
    smoothmin(d.x, da, .7, d.x);
    
    lfnoise_edge(34.*x.xy+3337.1, n.x);
    lfnoise_edge(14.*x.xy+4337.2, n.y);
    lfnoise_edge(54.*x.xy+5337.3, n.z);
    
    stroke(-2.*n.x+2.*n.y+n.z, .02, da);
    
    zextrude(x.z, -da, .02*n.x, da);
    stroke(da,mix(.2,.8,iScale),da);
    
    d = mix(d, vec2(da, 2.), step(da,d.x));
    
    d.x = max(d.x, x.z-.1+.02*n.z);
    
    vec2 ind,ia;
    dvoronoi(12.*x.xy, da, ind);
    vec2 y = x.xy-ind/12.;
    float r;
    rand(ind, r);
    r = length(y)-.02*r;
    
    d = mix(d,vec2(r,3.), step(r,d.x));
}

// Normal
const float dx = 5.e-4;
void normal(in vec3 x, out vec3 n);

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
    
    uv *= mix(1.,.1,iDial7);
    
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
            d += min(s.x,.005);
        }

        // Illumination
        l = normalize(x+c.yxx);
        if(i<N)
        {
            normal(x,n);
            
            vec3 ground = vec3(0.58,0.13,0.91);
            
            float na, nal;
			rand(floor(.33*iTime)*c.xx, na);
            rand(floor(.33*iTime)*c.xx+1., nal);
            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));
            
            float da;
            lfnoise_edge(21.*(x.xy+.1*mix(1.,5.,iDial0)*iTime*c.yx), da);
            stroke(da, .4, da);
            stroke(da-.1,.2,da);
            ground = mix(ground, 1.6*ground, step(0.,da));
            stroke(da-.01,.06,da);
            
            ground = mix(ground, .3*ground, step(0.,da));
            
            if(s.y == 1.)
            {
            	col = mix(ground, vec3(0.82,0.71,0.50), step(.002,x.z));
            	col = mix(col, 1.6*vec3(0.14,0.14,0.12), step(.018,x.z));
            }
            else if(s.y == 2.)
            {
                col = mix(ground, 1.6*vec3(0.91,0.84,0.79), step(.001,x.z));
            	col = mix(col, 1.6*vec3(0.95,0.80,0.68), step(.018,x.z));
            }
			else if(s.y == 3.)
            {
                col = mix(ground, vec3(0.71,0.00,0.17), step(.001,x.z));
            	col = mix(col, 1.6*vec3(0.88,0.77,0.99), step(.018,x.z));
            }
            mat3 RR;
            rot3(na*1.e3*vec3(1.1,1.5,1.9)+13.*length(col),RR);

            col = mix(col,.3*abs(RR*col),.5+.5*sin(1.*length(x.xy)+length(col)));
            col = mix(col,abs(RR*col),.5+.5*cos(21.*x.z+length(col)));
        }
    }
    
    // Colorize
    col = .2*col
        + 1.3*col*abs(dot(l,n))
        +.4*col*abs(pow(dot(reflect(-l,n),dir),3.));
    
    float dd;
    rand(1200.*uv, dd);
    col += dd*.1*c.xxx;
    
    col = mix(col, length(col)*c.xxx/sqrt(3.), .5);
    
    fragColor = clamp(vec4(col,1.0),0.,1.);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
