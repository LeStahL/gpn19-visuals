/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
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
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#version 130

uniform float iTime;
uniform float iFFTWidth;
uniform float iScale;
uniform float iHighScale;
uniform float iNBeats;
uniform vec2 iResolution;
uniform sampler1D iFFT;

out vec4 gl_FragColor;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

void rand(in vec2 x, out float num);
void lfnoise(in vec2 t, out float num);
void mfnoise(in vec2 x, in float fmin, in float fmax, in float alpha, out float num);

void camerasetup(in vec2 uv, out vec3 ro, out vec3 dir)
{
    vec3 right = c.xyy, up = c.yxy, target = c.yyy;
    ro = c.yyx+.3*vec3(cos(iTime), sin(iTime), 0.)*(1.-smoothstep(11., 13., iTime));
    dir = normalize(target + uv.x * right + uv.y * up - ro);
}

void stroke(in float d0, in float s, out float d);
void dcirclesegment(in vec2 x, in float R, in float p0, in float p1, out float d);
void dcircle(in vec2 x, in float R, out float d);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
void dlogo210(in vec2 x, in float R, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
void dbox3(in vec3 x, in vec3 b, out float d);
void dhexagonpattern(in vec2 p, out float d, out vec2 ind);
void rot3(in vec3 p, out mat3 rot);

// graph traversal for 210 logo effect
void textpre(in vec3 x, out vec2 sdf)
{
    float blend = 1.;
    dbox3(x, vec3(1., .5, .01+blend), sdf.x);
}

// Perform raymarching for bounding object
void marchbounds(in vec3 ro, in vec3 dir, in int N, in float eps, out vec3 x, out vec2 s, out float d, out bool flag)
{
    flag = false;
    for(int ia=0; ia<max(N,0); ++ia)
	{
        x = ro + d*dir;
        textpre(x,s);
        if(s.x < eps)
        {
            flag = true;
            break;
        }
        d += s.x;
	}
}

// 3D Effect on text in intro (210 logo)
void texteffect(in vec3 x, out vec2 sdf)
{
    // Start with z=0 plane
    sdf = vec2(x.z, 7.0);
    vec2 ind;
    float hex;
    dhexagonpattern(48.0*x.xy, hex, ind);
    
    // compute hexagon indices in cartesian coordinates
    vec2 cind = ind/48.0;
    
    // build up team210 logo (t < 12.)
    float inner_logo, logo_border; 
    dlogo210(3.5*cind, 1., inner_logo);
    stroke(inner_logo, 0.35, inner_logo);

    vec2 n;
    lfnoise(ind-2.*iTime, n.x);
    lfnoise(ind-2.*iTime-1337., n.y);
    
    // blend back to structure (t < 16., t > 12.)
    float blend = clamp(.5+.5*n.x, 0., 1.);
    
    if(inner_logo < 0.0 && blend >= 1.0e-3)
    {
        float noise;
        lfnoise(24.0*cind.xy-iTime, noise);
        zextrude(x.z,
                 1.5*x.z-inner_logo, 
                 .5*(0.5+0.5*noise)*blend*(.5+.5*iScale),
                 sdf.x);
        stroke(sdf.x, 0.05*blend, sdf.x);
        sdf.y = 7.0;
    }
    stroke(sdf.x,0.1,sdf.x);
    
//     // Add guard objects for debugging
//     float dr = .02;
//     vec3 y = mod(x,dr)-.5*dr;
//     float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
//     guard = abs(guard)+dr*.1;
//     sdf.x = min(sdf.x, guard);
}

// Perform raymarching for actual object
void marchscene(in vec3 ro, in vec3 dir, in int N, in float eps, out vec3 x, out vec2 s, out float d, out bool flag)
{
    flag = false;
    for(int ia=0; ia<max(N,0); ++ia)
	{
        x = ro + d*dir;
        texteffect(x,s);
        if(s.x < eps)
        {
            flag = true;
            break;
        }
        d += min(s.x,.005);
	}
}

void calcnormal(in vec3 x, in float eps, out vec3 n)
{
    vec2 s, sp;
    texteffect(x, s);
    texteffect(x+eps*c.xyy, sp);
    n.x = sp.x-s.x;
    texteffect(x+eps*c.yxy, sp);
    n.y = sp.x-s.x;
    texteffect(x+eps*c.yyx, sp);
    n.z = sp.x-s.x;
    n = normalize(n);
}

// Initial intro
void background2(in vec2 uv, out vec3 col)
{
    col = c.yyy;
    
    // hexagonal grid
    float d, d0;
    vec2 ind;
    dhexagonpattern(48.0*uv, d0, ind);
    d = -d0;
    stroke(d, 0.1, d);
    vec2 cind = ind/48.0;
    
    // build up team210 logo (t < 12.)
    float inner_logo, logo_border; 
    dlogo210(3.5*cind, 1., inner_logo);
    stroke(inner_logo, 0.35, inner_logo);
    stroke(inner_logo, 0.08, logo_border);
    
    vec2 n;
    lfnoise(ind-2.*iTime, n.x);
    lfnoise(ind-2.*iTime-1337., n.y);
    
    // blend back to structure (t < 16., t > 12.)
    float blend = 0.;
    
    inner_logo = mix(inner_logo, d0, blend);
    logo_border = mix(logo_border, d0, blend);

    // make background change the color with time
    vec2 dt;
    lfnoise(15.0*cind+2.0, dt.x);
    lfnoise(15.0*cind+3.0, dt.y);
    dt *= 2.;
    float dm, dm2;
    lfnoise(50.0*cind, dm);
    dm = 0.5+0.5*dm;
    lfnoise(6.5*cind-dt-2.0*iTime*c.xx, dm2);
    dm2 = 0.5+0.5*dm2;
    
    // Colors
    vec3 orange = vec3(0.26,0.35,0.40);
    orange = mix(vec3(.60,0.24,0.60), orange, dm2);
    vec3 gray = .25*length(orange)/sqrt(3.)*c.xxx;
    
    float mat = n.x;
//     rand(i*c.xx, mat);
    float phi = atan(uv.y, uv.x),
        dhex,
        na,
        nal;
    rand(floor(.33*iTime)*c.xx, na);
    rand(floor(.33*iTime)*c.xx+1., nal);
    na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));
    
    mat3 RRR;
    rot3(na*1.e3*vec3(1.1,1.5,1.9),RRR);

    orange = mix((.5+.5*n.y)*orange,(1.+.8*mat)*abs(RRR*orange),.5+.5*n.x);
    orange = mix(orange,.2*orange,.5+.5*n.y);
  
    col = orange;
    col = mix(col, gray, step(0.,inner_logo));
    
    col = mix(col, mix(gray,orange,step(abs(d0)-.4,0.)), step(-abs(d0)+.2,0.));
    col = mix(col, 4.*orange, step(logo_border,0.));
    
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), s = c.xy;

    if(length(uv) > .5)
    {
        fragColor = vec4(c.yyy, 0.);
        return;
    }
    uv *= 1.7;
    
	vec3 ro, x, dir;
    
    float d = 0.;
    bool hit = false;
    
    vec3 col = c.yyy;
                
	camerasetup(uv, ro, dir);
    d = (.5-ro.z)/dir.z;
    marchbounds(ro, dir, 150, 2.0e-4, x, s, d, hit);

    if(hit) hit = false;
    else d = -ro.z/dir.z;
    marchscene(ro, dir, 800, 2.0e-4, x, s, d, hit);
    
    if(hit)
    {
        vec3 n;
        calcnormal(x, 2.0e-4, n);

        float rs = 1.9;
        vec3 l = x+1.*c.yyx,
        	re = normalize(reflect(-l,n));
        float rev = abs(dot(re,dir)), ln = abs(dot(l,n));
		background2(x.xy, col);
        col = mix(col, .5*col, clamp(x.z/.2,0.,1.));
    }
    else
        background2((ro-ro.z/dir.z*dir).xy, col);

    col = clamp(col, 0., 1.);
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
