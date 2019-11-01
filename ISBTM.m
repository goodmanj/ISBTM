function out = model
%
% ISBTM.m
%
% Model exported on Nov 1 2019, 06:05 by COMSOL 5.4.0.388.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\goodm\OneDrive\Documents\ISBTM');

model.label('ISBTM.mph');

model.param.set('shellwidth', '3e4[m]');
model.param.set('shellthick', '1e4[m]');
model.param.set('bumpwidth', '5000[m]');
model.param.set('bumpheight', '-2500[m]');
model.param.group.create('par3');
model.param('par3').set('g', '1.3 [m/s^2]', 'Europa gravity');
model.param.group.create('par2');
model.param('par2').set('dd', '1e-3 [m]', 'Grain size meters');
model.param('par2').set('minsr', '1e-16[1/s]', 'Minimum shear rate magnitude for viscosity calculations');
model.param('par2').set('GlenA1', '3.61e-13 [Pa^(-3)/s]', 'Macayeal, page 399');
model.param('par2').set('GlenA2', '1.73e3 [Pa^(-3)/s]', 'Macayeal page 399');
model.param('par2').set('GlenQ1', '60000 [J/mol]', 'Macayeal 399');
model.param('par2').set('GlenQ2', '139000 [J/mol]', 'Macayeal pg 399');
model.param('par2').set('A1', '1.2e-10', 'diffusion const, si units');
model.param('par2').set('A2', '2.2e-7', 'basal const, si units');
model.param('par2').set('A3', '6.2e-14', 'gbs const, si units');
model.param('par2').set('A4', '4e-19', 'disl constant, si units');
model.param('par2').set('n1', '1', 'diff exponent');
model.param('par2').set('n2', '2.4', 'basal exponent');
model.param('par2').set('n3', '1.8', 'gbs exponent');
model.param('par2').set('n4', '4', 'disl exponent');
model.param('par2').set('p1', '2', 'diff grain-size dependence');
model.param('par2').set('p2', '0', 'basal grain-size dependence');
model.param('par2').set('p3', '1.4', 'gbs grain-size dependence');
model.param('par2').set('p4', '0', 'disl grain-size dependence');
model.param('par2').set('Q1', '59400 [J/mol]', 'diff activation energy');
model.param('par2').set('Q2', '60000[J/mol]', 'basal activation energy');
model.param('par2').set('Q3', '49000[J/mol]', 'gbs activation energy');
model.param('par2').set('Q4', '60000[J/mol]', 'disl activation energy');
model.param('par2').set('RR', '8.31 [J/(mol*K)]', 'ideal gas const');
model.param('par2').set('A3w', '1.11e16', 'enhanced GBS warm ice');
model.param('par2').set('A4w', '8e5', 'TWEAKED enhanced dislocation warm');
model.param('par2').set('Q3w', '192000 [J/mol]', 'enhanced gbs Q');
model.param('par2').set('Q4w', '180000 [J/mol]', 'enhanced disl Q');
model.param('par2').set('Tgbs', '255 [K]', 'threshhold temp for enhanced gbs');
model.param('par2').set('Tdisl', '258 [K]', 'threshhold temp for enhanced disl');
model.param.group.create('par4');
model.param('par4').set('rho2', '1000 [kg/m^3]', 'Density of underlying fluid');
model.param('par4').set('mu2', '1e10 [Pa*s]', 'Viscosity of underlying fluid');
model.param.group.create('par5');
model.param('par5').set('Qocean', '.065 [W/m^2]', 'Average heat flux from ocean');
model.param('par5').set('Tbottom', '273.15 [K]', 'Temperature of Bottom of ice cap');
model.param('par5').set('Tstart', '200[K]', 'Starting temperature of ice block');
model.param('par5').set('Ttop', '100[K]', 'Temperature of Top of Europa (Vacuum/Ice Barrier)');
model.param.group.create('par6');
model.param('par6').set('Cal_per_Mol_K_to_J_per_Kg_K', '4184/18.015');
model.param('par6').set('L', '333.4e3 [J/kg]', 'Latent Heat of Fusion for Water');
model.param.label('Geometry Parameters');
model.param('par3').label('Planetary Parameters');
model.param('par2').label('Flow law parameters');
model.param('par2').comments(['Goldsby flow law parameters: Warm dislocation creep constant increased from 6e4 (Goldsby 2001 LPSC abstract value) to 8e5 to ensure continuity with cold value, and better match Goldsby 2001 figure.\n']);
model.param('par4').label('Sticky Water Parameters');
model.param('par5').label('Heat Flow Parameters');
model.param('par6').label('Ice Thermodynamic Parameters');
model.param('par6').comments(['CRC Handbook of Physics and Chemistry, 97th edition, page 6-12\n']);

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.func.create('an1', 'Analytic');
model.func.create('an2', 'Analytic');
model.func.create('int1', 'Interpolation');
model.func.create('int2', 'Interpolation');
model.func.create('int3', 'Interpolation');
model.func.create('an3', 'Analytic');
model.func.create('an4', 'Analytic');
model.func.create('an5', 'Analytic');
model.func.create('an8', 'Analytic');
model.func.create('an6', 'Analytic');
model.func.create('an9', 'Analytic');
model.func.create('an7', 'Analytic');
model.func('an1').label('Glen Flow Law A Parameter');
model.func('an1').set('funcname', 'GlenA');
model.func('an1').set('expr', 'GlenA1*exp(-GlenQ1/(R_const*T))*(T <= 263.15) + GlenA2*exp(-GlenQ2/(R_const*T))*(T >= 263.15)');
model.func('an1').set('args', {'T'});
model.func('an1').set('argunit', 'K');
model.func('an1').set('plotargs', {'T' '100' '273'});
model.func('an2').label('Glen Viscosity');
model.func('an2').set('funcname', 'mu_glen');
model.func('an2').set('expr', '(GlenA(T )* 8 * (sr)^2 )^-(1/3)');
model.func('an2').set('args', {'T' 'sr'});
model.func('an2').set('argunit', 'K,1/s');
model.func('an2').set('fununit', 'Pa*s');
model.func('an2').set('plotargs', {'T' '100' '273'; 'sr' '1e-3' '1e-8'});
model.func('int1').label('Ice Specific Heat');
model.func('int1').comments(['CRC Handbook of Physics and Chemistry, 97th edition, page 6-12\n']);
model.func('int1').set('funcname', 'Cpice');
model.func('int1').set('table', {'13' '36';  ...
'33' '270';  ...
'53' '470';  ...
'73' '650';  ...
'93' '820';  ...
'113' '970';  ...
'133' '1110';  ...
'153' '1250';  ...
'173' '1380';  ...
'193' '1520';  ...
'213' '1660';  ...
'223' '1730';  ...
'233' '1800';  ...
'243' '1880';  ...
'253' '1950';  ...
'263' '2020';  ...
'273' '2100'});
model.func('int1').set('interp', 'piecewisecubic');
model.func('int1').set('extrap', 'linear');
model.func('int1').set('argunit', 'K');
model.func('int1').set('fununit', 'J/(kg*K)');
model.func('int2').label('Ice Density');
model.func('int2').comments(['CRC Handbook of Physics and Chemistry, 97th edition, page 6-12\n']);
model.func('int2').set('funcname', 'rhoice');
model.func('int2').set('table', {'13' '933.8';  ...
'33' '933.8';  ...
'53' '933.7';  ...
'73' '933.6';  ...
'93' '933.2';  ...
'113' '932.6';  ...
'133' '931.7';  ...
'153' '930.4';  ...
'173' '928.8';  ...
'193' '926.9';  ...
'213' '924.7';  ...
'223' '923.5';  ...
'233' '922.2';  ...
'243' '920.9';  ...
'253' '919.6';  ...
'263' '918.2';  ...
'273' '916.7'});
model.func('int2').set('interp', 'piecewisecubic');
model.func('int2').set('extrap', 'linear');
model.func('int2').set('argunit', 'K');
model.func('int2').set('fununit', 'kg/m^3');
model.func('int3').label('Ice Thermal Conductivity');
model.func('int3').comments(['CRC Handbook of Physics and Chemistry, 97th edition, page 6-12\n']);
model.func('int3').set('funcname', 'lambdaice');
model.func('int3').set('table', {'33' '20';  ...
'53' '12.2';  ...
'73' '8.9';  ...
'93' '7.0';  ...
'113' '5.7';  ...
'133' '4.9';  ...
'153' '4.2';  ...
'173' '3.69';  ...
'193' '3.27';  ...
'213' '2.93';  ...
'223' '2.77';  ...
'233' '2.63';  ...
'243' '2.50';  ...
'253' '2.38';  ...
'263' '2.26';  ...
'273' '2.16'});
model.func('int3').set('interp', 'piecewisecubic');
model.func('int3').set('extrap', 'linear');
model.func('int3').set('argunit', 'K');
model.func('int3').set('fununit', 'W/(m*K)');
model.func('an3').label('Golds_diffusion');
model.func('an3').set('funcname', 'eta1');
model.func('an3').set('expr', '(dd^(p1/n1))/((2*A1)^(1/n1))*exp(Q1/(RR*T*n1))*sr^((1-n1)/n1)');
model.func('an3').set('args', {'T' 'sr'});
model.func('an3').set('argunit', 'K,1/s');
model.func('an3').set('fununit', 'Pa*s');
model.func('an3').set('plotargs', {'T' '100' '273'; 'sr' '1e-10' '1e-8'});
model.func('an4').label('Golds_basal');
model.func('an4').set('funcname', 'eta2');
model.func('an4').set('expr', '(dd^(p2/n2))/((2*A2)^(1/n2))*exp(Q2/(RR*T*n2))*sr^((1-n2)/n2)');
model.func('an4').set('args', {'T' 'sr'});
model.func('an4').set('argunit', 'K,1/s');
model.func('an4').set('fununit', 'Pa*s');
model.func('an4').set('plotargs', {'T' '100' '273'; 'sr' '1e-10' '1e-8'});
model.func('an5').label('Golds_gbs');
model.func('an5').set('funcname', 'eta3');
model.func('an5').set('expr', '(dd^(p3/n3))/((2*A3)^(1/n3))*exp(Q3/(RR*T*n3))*sr^((1-n3)/n3)');
model.func('an5').set('args', {'T' 'sr'});
model.func('an5').set('argunit', 'K,1/s');
model.func('an5').set('fununit', 'Pa*s');
model.func('an5').set('plotargs', {'T' '100' '273'; 'sr' '1e-10' '1e-8'});
model.func('an8').label('Golds_gbs warm');
model.func('an8').set('funcname', 'eta3w');
model.func('an8').set('expr', '(T<=Tgbs)*((dd^(p3/n3))/((2*A3)^(1/n3)))*exp(Q3/(RR*T*n3))*sr^((1-n3)/n3)+(T>Tgbs)*(dd^(p3/n3))/((2*A3w)^(1/n3))*exp(Q3w/(RR*T*n3))*sr^((1-n3)/n3)');
model.func('an8').set('args', {'T' 'sr'});
model.func('an8').set('argunit', 'K,1/s');
model.func('an8').set('fununit', 'Pa*s');
model.func('an8').set('plotargs', {'T' '100' '273'; 'sr' '1e-10' '1e-8'});
model.func('an6').label('Golds_dislocation');
model.func('an6').set('funcname', 'eta4');
model.func('an6').set('expr', '(dd^(p4/n4))/((2*A4)^(1/n4))*exp(Q4/(RR*T*n4))*sr^((1-n4)/n4)');
model.func('an6').set('args', {'T' 'sr'});
model.func('an6').set('argunit', 'K,1/s');
model.func('an6').set('fununit', 'Pa*s');
model.func('an6').set('plotargs', {'T' '100' '273'; 'sr' '1e-10' '1e-8'});
model.func('an9').label('Golds_dislocation warm');
model.func('an9').set('funcname', 'eta4w');
model.func('an9').set('expr', '(T<=Tdisl)*(1)/((2*A4)^(1/n4))*exp(Q4/(RR*T*n4))*sr^((1-n4)/n4)+(T>Tdisl)*(1)/((2*A4w)^(1/n4))*exp(Q4w/(RR*T*n4))*sr^((1-n4)/n4)');
model.func('an9').set('args', {'T' 'sr'});
model.func('an9').set('argunit', 'K,1/s');
model.func('an9').set('fununit', 'Pa*s');
model.func('an9').set('plotargs', {'T' '100' '273'; 'sr' '1e-10' '1e-8'});
model.func('an7').label('Goldsby viscosity');
model.func('an7').set('funcname', 'mu_goldsby');
model.func('an7').set('expr', 'min(min(eta1(T,sr),max(eta2(T,sr),eta3w(T,sr))),eta4w(T,sr))');
model.func('an7').set('args', {'T' 'sr'});
model.func('an7').set('argunit', 'K, 1/s');
model.func('an7').set('fununit', 'Pa*s');
model.func('an7').set('plotargs', {'T' '100' '270'; 'sr' '1e-10' '1e-8'});

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('pc1', 'ParametricCurve');
model.component('comp1').geom('geom1').feature('pc1').set('pos', {'0' '-shellthick-bumpheight*exp(-(shellwidth/(2*bumpwidth))^2)'});
model.component('comp1').geom('geom1').feature('pc1').set('parmin', '-shellwidth/2');
model.component('comp1').geom('geom1').feature('pc1').set('parmax', 'shellwidth/2');
model.component('comp1').geom('geom1').feature('pc1').set('coord', {'s' 'bumpheight*exp(-(s[m]/bumpwidth)^2)'});
model.component('comp1').geom('geom1').create('b1', 'BezierPolygon');
model.component('comp1').geom('geom1').feature('b1').set('p', {'-shellwidth/2' '-shellwidth/2'; '-shellthick' '0'});
model.component('comp1').geom('geom1').feature('b1').set('degree', [1]);
model.component('comp1').geom('geom1').feature('b1').set('w', [1 1]);
model.component('comp1').geom('geom1').create('b2', 'BezierPolygon');
model.component('comp1').geom('geom1').feature('b2').set('p', {'-shellwidth/2' 'shellwidth/2'; '0' '0'});
model.component('comp1').geom('geom1').feature('b2').set('degree', [1]);
model.component('comp1').geom('geom1').feature('b2').set('w', [1 1]);
model.component('comp1').geom('geom1').create('b3', 'BezierPolygon');
model.component('comp1').geom('geom1').feature('b3').set('p', {'shellwidth/2' 'shellwidth/2'; '0' '-shellthick'});
model.component('comp1').geom('geom1').feature('b3').set('degree', [1]);
model.component('comp1').geom('geom1').feature('b3').set('w', [1 1]);
model.component('comp1').geom('geom1').create('csol1', 'ConvertToSolid');
model.component('comp1').geom('geom1').feature('csol1').selection('input').set({'b1' 'b2' 'b3' 'pc1'});
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'-shellwidth/2' '-(shellthick*1.2+abs(bumpheight))'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'shellwidth' 'shellthick*1.2+abs(bumpheight)'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.view.create('view2', 3);
model.view.create('view3', 3);

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat1').selection.set([2]);
model.component('comp1').material('mat2').selection.set([1]);

model.component('comp1').common.create('free1', 'DeformingDomain');
model.component('comp1').common.create('pnmv1', 'PrescribedNormalMeshVelocity');
model.component('comp1').common.create('pnmv2', 'PrescribedNormalMeshVelocity');
model.component('comp1').common.create('slip1', 'Slip');
model.component('comp1').common('free1').selection.set([1 2]);
model.component('comp1').common('pnmv1').selection.set([4]);
model.component('comp1').common('pnmv2').selection.set([7]);
model.component('comp1').common('slip1').selection.set([1 3 5 6]);

model.component('comp1').physics.create('spf', 'LaminarFlow', 'geom1');
model.component('comp1').physics('spf').create('vf1', 'VolumeForce', 2);
model.component('comp1').physics('spf').feature('vf1').selection.set([1 2]);
model.component('comp1').physics('spf').create('bs2', 'BoundaryStress', 1);
model.component('comp1').physics('spf').feature('bs2').selection.set([4]);
model.component('comp1').physics('spf').create('prpc1', 'PressurePointConstraint', 0);
model.component('comp1').physics('spf').feature('prpc1').selection.set([3]);
model.component('comp1').physics.create('ht', 'HeatTransferInFluids', 'geom1');
model.component('comp1').physics('ht').selection.set([2]);
model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp1').selection.set([4]);
model.component('comp1').physics('ht').create('temp2', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp2').selection.set([7]);

model.component('comp1').multiphysics.create('fc1', 'FlowCoupling', -1);
model.component('comp1').multiphysics.create('tc1', 'TemperatureCoupling', -1);

model.component('comp1').mesh('mesh1').create('bl1', 'BndLayer');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('bl1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('bl1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('bl1').create('blp', 'BndLayerProp');
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').selection.set([7]);

model.component('comp1').view('view1').axis.set('xmin', -18840);
model.component('comp1').view('view1').axis.set('xmax', 18840);
model.component('comp1').view('view1').axis.set('ymin', -22219.998046875);
model.component('comp1').view('view1').axis.set('ymax', 7720.00048828125);
model.component('comp1').view('view1').axis.set('xextra', 0);
model.component('comp1').view('view1').axis.set('yextra', 0);

model.component('comp1').material('mat1').label('Ice');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rhoice(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'mu_goldsby(T,spf.sr)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'lambdaice(T)' '0' '0' '0' 'lambdaice(T)' '0' '0' '0' 'lambdaice(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cpice(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1');
model.component('comp1').material('mat2').label('Sticky Water');
model.component('comp1').material('mat2').propertyGroup('def').set('density', 'rho2');
model.component('comp1').material('mat2').propertyGroup('def').set('dynamicviscosity', 'mu2');

model.component('comp1').common('free1').set('smoothingType', 'hyperelastic');
model.component('comp1').common('pnmv1').label('Upper Normal Mesh Velocity');
model.component('comp1').common('pnmv1').set('prescribedNormalVelocity', 'u*nx+v*ny');
model.component('comp1').common('pnmv2').label('Lower Normal Mesh Velocity 1');
model.component('comp1').common('pnmv2').comments(['Bottom surface of ice moves in response to both ice flow and melting/freezing.  Melt velocity equals excess heat delivered from ocean divided by latent heat of fusion (in J/meter thickness)\n']);
model.component('comp1').common('pnmv2').set('prescribedNormalVelocity', 'u*nx+(v+(Qocean-ht.dfluxy)/(L*ht.rho))*ny');

model.component('comp1').physics('spf').feature('fp1').featureInfo('info').set('spf.sr', {'sqrt(2*spf.srijxx^2+2*spf.srijxy^2+2*spf.srijxz^2+2*spf.srijyx^2+2*spf.srijyy^2+2*spf.srijyz^2+2*spf.srijzx^2+2*spf.srijzy^2+2*spf.srijzz^2+minsr^2)'});
model.component('comp1').physics('spf').feature('wallbc1').set('BoundaryCondition', 'Slip');
model.component('comp1').physics('spf').feature('vf1').set('F', {'0'; '-spf.rho*g'; '0'});
model.component('comp1').physics('spf').feature('prpc1').set('p0', 1);
model.component('comp1').physics('ht').feature('temp1').set('T0', 100);
model.component('comp1').physics('ht').feature('temp2').set('T0', '273.15[K]');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 4);
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('bl1').active(false);
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').set('blnlayers', 3);
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').set('inittype', 'blhmin');
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').set('blhmin', 250);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('stat').set('activate', {'spf' 'off' 'ht' 'on' 'frame:spatial1' 'off'});
model.study('std1').feature('stat').set('activateCoupling', {'fc1' 'off' 'tc1' 'off'});

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').create('su1', 'StoreSolution');
model.sol('sol1').create('st2', 'StudyStep');
model.sol('sol1').create('v2', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').create('i2', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').feature('t1').create('se1', 'Segregated');
model.sol('sol1').feature('t1').create('d1', 'Direct');
model.sol('sol1').feature('t1').create('d2', 'Direct');
model.sol('sol1').feature('t1').create('d3', 'Direct');
model.sol('sol1').feature('t1').create('i1', 'Iterative');
model.sol('sol1').feature('t1').create('i2', 'Iterative');
model.sol('sol1').feature('t1').feature('se1').create('ss1', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').create('ss2', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').create('ss3', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').create('ll1', 'LowerLimit');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ssDef');
model.sol('sol1').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');

model.result.dataset.create('edg1', 'Edge2D');
model.result.dataset.create('edg2', 'Edge2D');
model.result.dataset.create('cpt1', 'CutPoint2D');
model.result.dataset.create('cpt2', 'CutPoint2D');
model.result.dataset.create('grid1', 'Grid1D');
model.result.dataset.create('grid2', 'Grid1D');
model.result.dataset.create('an7_ds1', 'Grid2D');
model.result.dataset('edg1').selection.set([4]);
model.result.dataset('edg2').selection.set([7]);
model.result.dataset('grid1').set('data', 'none');
model.result.dataset('grid2').set('data', 'none');
model.result.dataset('an7_ds1').set('data', 'none');
model.result.dataset.remove('dset2');
model.result.dataset.remove('dset3');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg8', 'PlotGroup2D');
model.result.create('pg10', 'PlotGroup2D');
model.result.create('pg11', 'PlotGroup2D');
model.result.create('pg15', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup1D');
model.result.create('pg4', 'PlotGroup1D');
model.result.create('pg7', 'PlotGroup1D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').create('mesh1', 'Mesh');
model.result('pg1').create('arws1', 'ArrowSurface');
model.result('pg2').create('con1', 'Contour');
model.result('pg8').create('surf1', 'Surface');
model.result('pg10').create('surf1', 'Surface');
model.result('pg11').create('surf1', 'Surface');
model.result('pg15').create('surf1', 'Surface');
model.result('pg3').create('lngr1', 'LineGraph');
model.result('pg4').create('lngr1', 'LineGraph');
model.result('pg7').create('ptgr1', 'PointGraph');
model.result('pg7').create('ptgr2', 'PointGraph');
model.result.export.create('anim1', 'Animation');

model.study('std1').feature('time').set('tunit', 'a');
model.study('std1').feature('time').set('tlist', 'range(0,100,20000)');
model.study('std1').feature('time').set('plot', true);
model.study('std1').feature('time').set('plotfreq', 'tsteps');
model.study('std1').feature('time').set('useinitsol', true);
model.study('std1').feature('time').set('initmethod', 'sol');
model.study('std1').feature('time').set('initstudy', 'std1');
model.study('std1').feature('time').set('initsoluse', 'sol8');
model.study('std1').feature('time').set('solnum', 'auto');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {'20.0[a]'});
model.sol('sol1').feature('v1').feature('comp1_spatial_disp').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_spatial_disp').set('scaleval', 424.9712931481372);
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('s1').feature('fc1').set('termonres', false);
model.sol('sol1').feature('s1').feature('d1').label('PARDISO (ht)');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s1').feature('i1').label('Algebraic Multigrid (ht)');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i2').label('Geometric Multigrid (ht)');
model.sol('sol1').feature('s1').feature('i2').set('rhob', 20);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('st2').set('studystep', 'time');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('initsoluse', 'sol8');
model.sol('sol1').feature('v2').set('solnum', 'auto');
model.sol('sol1').feature('v2').set('resscalemethod', 'manual');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('notsolnum', 'auto');
model.sol('sol1').feature('v2').set('clist', {'range(0,100,20000)' '20.0[a]'});
model.sol('sol1').feature('v2').feature('comp1_spatial_disp').set('scalemethod', 'manual');
model.sol('sol1').feature('v2').feature('comp1_spatial_disp').set('scaleval', 424.9712931481372);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').set('tunit', 'a');
model.sol('sol1').feature('t1').set('tlist', 'range(0,100,20000)');
model.sol('sol1').feature('t1').set('rtol', 0.005);
model.sol('sol1').feature('t1').set('atolglobalfactor', 0.05);
model.sol('sol1').feature('t1').set('atolmethod', {'comp1_p' 'scaled' 'comp1_spatial_disp' 'global' 'comp1_spatial_lm' 'global' 'comp1_spatial_lm_nv' 'global' 'comp1_T' 'global'  ...
'comp1_u' 'global'});
model.sol('sol1').feature('t1').set('atolfactor', {'comp1_p' '1' 'comp1_spatial_disp' '0.1' 'comp1_spatial_lm' '0.1' 'comp1_spatial_lm_nv' '0.1' 'comp1_T' '0.1'  ...
'comp1_u' '0.1'});
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('stabcntrl', true);
model.sol('sol1').feature('t1').set('bwinitstepfrac', 0.01);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('plot', true);
model.sol('sol1').feature('t1').set('plotfreq', 'tsteps');
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').feature('se1').set('ntolfact', 0.5);
model.sol('sol1').feature('t1').feature('se1').set('segstabacc', 'segaacc');
model.sol('sol1').feature('t1').feature('se1').set('segaaccdim', 5);
model.sol('sol1').feature('t1').feature('se1').set('segaaccmix', 0.9);
model.sol('sol1').feature('t1').feature('se1').set('segaaccdelay', 1);
model.sol('sol1').feature('t1').feature('se1').feature('ss1').label('Heat transfer T');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('linsolver', 'd1');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('subdamp', 0.7);
model.sol('sol1').feature('t1').feature('se1').feature('ss2').label('Velocity u, Pressure p');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('segvar', {'comp1_u' 'comp1_p'});
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('linsolver', 'd2');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('subdamp', 0.8);
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('subjtech', 'once');
model.sol('sol1').feature('t1').feature('se1').feature('ss3').label('Spatial mesh displacement');
model.sol('sol1').feature('t1').feature('se1').feature('ss3').set('segvar', {'comp1_spatial_disp' 'comp1_spatial_lm' 'comp1_spatial_lm_nv'});
model.sol('sol1').feature('t1').feature('se1').feature('ss3').set('linsolver', 'd3');
model.sol('sol1').feature('t1').feature('se1').feature('ss3').set('subdamp', 0.8);
model.sol('sol1').feature('t1').feature('se1').feature('ss3').set('subjtech', 'once');
model.sol('sol1').feature('t1').feature('se1').feature('ll1').set('lowerlimit', 'comp1.T 0');
model.sol('sol1').feature('t1').feature('d2').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d2').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('d3').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d3').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('i1').label('Algebraic Multigrid (ht)');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('i2').label('Geometric Multigrid (ht)');
model.sol('sol1').feature('t1').feature('i2').set('rhob', 20);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

model.result.dataset('edg1').label('Top surface');
model.result.dataset('edg2').label('Bottom');
model.result.dataset('cpt1').label('Top center');
model.result.dataset('cpt1').set('pointx', 0);
model.result.dataset('cpt1').set('pointy', 0);
model.result.dataset('cpt1').set('bndsnap', true);
model.result.dataset('cpt2').label('Bottom center');
model.result.dataset('cpt2').set('pointx', 0);
model.result.dataset('cpt2').set('pointy', '-shellthick+bumpheight');
model.result.dataset('cpt2').set('bndsnap', true);
model.result.dataset('cpt2').set('pointvar', 'cpt1n');
model.result.dataset('grid1').label('Temperature range for gbs');
model.result.dataset('grid1').set('function', 'an8');
model.result.dataset('grid1').set('par1', 'T');
model.result.dataset('grid1').set('parmin1', 100);
model.result.dataset('grid1').set('parmax1', 273);
model.result.dataset('grid2').label('Temperature range for dislocation');
model.result.dataset('grid2').set('function', 'an9');
model.result.dataset('grid2').set('par1', 'T');
model.result.dataset('grid2').set('parmin1', 100);
model.result.dataset('grid2').set('parmax1', 273);
model.result.dataset('an7_ds1').set('function', 'all');
model.result.dataset('an7_ds1').set('par1', 'T');
model.result.dataset('an7_ds1').set('parmin1', 100);
model.result.dataset('an7_ds1').set('parmax1', 270);
model.result.dataset('an7_ds1').set('par2', 'sr');
model.result.dataset('an7_ds1').set('parmin2', 1.0E-10);
model.result.dataset('an7_ds1').set('parmax2', 1.0E-8);
model.result('pg1').label('Velocity (spf)');
model.result('pg1').set('looplevel', [197]);
model.result('pg1').set('view', 'view1');
model.result('pg1').set('edgecolor', 'white');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('expr', 'spf.U * (spf.rho < 950)');
model.result('pg1').feature('surf1').set('descr', 'spf.U * (spf.rho < 950)');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('mesh1').active(false);
model.result('pg1').feature('mesh1').set('elemcolor', 'none');
model.result('pg1').feature('arws1').set('scale', 2.83147998450761E11);
model.result('pg1').feature('arws1').set('color', 'white');
model.result('pg1').feature('arws1').set('scaleactive', false);
model.result('pg2').label('Pressure (spf)');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature('con1').label('Contour');
model.result('pg2').feature('con1').set('expr', 'p');
model.result('pg2').feature('con1').set('unit', 'Pa');
model.result('pg2').feature('con1').set('descr', 'Pressure');
model.result('pg2').feature('con1').set('number', 40);
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('resolution', 'normal');
model.result('pg8').label('Temperature');
model.result('pg8').set('frametype', 'spatial');
model.result('pg8').feature('surf1').set('expr', 'T');
model.result('pg8').feature('surf1').set('unit', 'K');
model.result('pg8').feature('surf1').set('descr', 'Temperature');
model.result('pg8').feature('surf1').set('resolution', 'normal');
model.result('pg10').label('Viscosity');
model.result('pg10').set('looplevel', [1]);
model.result('pg10').set('frametype', 'spatial');
model.result('pg10').feature('surf1').set('expr', 'log10(spf.mu)');
model.result('pg10').feature('surf1').set('unit', '');
model.result('pg10').feature('surf1').set('descr', 'log10(spf.mu)');
model.result('pg10').feature('surf1').set('resolution', 'normal');
model.result('pg11').label('Strain rate');
model.result('pg11').set('frametype', 'spatial');
model.result('pg11').feature('surf1').set('expr', 'log10(spf.sr)');
model.result('pg11').feature('surf1').set('unit', '');
model.result('pg11').feature('surf1').set('descr', 'log10(spf.sr)');
model.result('pg11').feature('surf1').set('resolution', 'normal');
model.result('pg15').label('Stress');
model.result('pg15').set('frametype', 'spatial');
model.result('pg15').feature('surf1').set('expr', 'log10(2*spf.mu*spf.sr)');
model.result('pg15').feature('surf1').set('unit', '');
model.result('pg15').feature('surf1').set('descr', 'log10(2*spf.mu*spf.sr)');
model.result('pg15').feature('surf1').set('resolution', 'normal');
model.result('pg3').label('Top Height');
model.result('pg3').set('data', 'edg1');
model.result('pg3').set('xlabel', 'Arc length (m)');
model.result('pg3').set('ylabel', 'y-coordinate (m)');
model.result('pg3').set('xlabelactive', false);
model.result('pg3').set('ylabelactive', false);
model.result('pg3').feature('lngr1').set('expr', 'y');
model.result('pg3').feature('lngr1').set('unit', 'm');
model.result('pg3').feature('lngr1').set('descr', 'y-coordinate');
model.result('pg3').feature('lngr1').set('titletype', 'custom');
model.result('pg3').feature('lngr1').set('expressionintitle', true);
model.result('pg3').feature('lngr1').set('resolution', 'normal');
model.result('pg4').label('Bottom Height');
model.result('pg4').set('data', 'edg2');
model.result('pg4').set('xlabel', 'Arc length (m)');
model.result('pg4').set('ylabel', 'y-coordinate (m)');
model.result('pg4').set('xlabelactive', false);
model.result('pg4').set('ylabelactive', false);
model.result('pg4').feature('lngr1').set('expr', 'y');
model.result('pg4').feature('lngr1').set('unit', 'm');
model.result('pg4').feature('lngr1').set('descr', 'y-coordinate');
model.result('pg4').feature('lngr1').set('resolution', 'normal');
model.result('pg7').label('Center height');
model.result('pg7').set('data', 'none');
model.result('pg7').set('xlabel', 'Time (a)');
model.result('pg7').set('xlabelactive', false);
model.result('pg7').feature('ptgr1').set('data', 'cpt1');
model.result('pg7').feature('ptgr1').set('looplevelinput', {'all'});
model.result('pg7').feature('ptgr1').set('expr', 'y');
model.result('pg7').feature('ptgr1').set('unit', 'm');
model.result('pg7').feature('ptgr1').set('descr', 'y-coordinate');
model.result('pg7').feature('ptgr1').set('legend', true);
model.result('pg7').feature('ptgr1').set('legendmethod', 'manual');
model.result('pg7').feature('ptgr1').set('legends', {'Top center height'});
model.result('pg7').feature('ptgr2').label('Height of bottom center');
model.result('pg7').feature('ptgr2').set('data', 'cpt2');
model.result('pg7').feature('ptgr2').set('looplevelinput', {'all'});
model.result('pg7').feature('ptgr2').set('expr', 'y+shellthick');
model.result('pg7').feature('ptgr2').set('unit', 'm');
model.result('pg7').feature('ptgr2').set('descr', 'y+shellthick');
model.result('pg7').feature('ptgr2').set('legend', true);
model.result('pg7').feature('ptgr2').set('legendmethod', 'manual');
model.result('pg7').feature('ptgr2').set('legends', {'Bottom center height' ''});
model.result.export('anim1').set('target', 'player');
model.result.export('anim1').set('showframe', 25);
model.result.export('anim1').set('shownparameter', '2');
model.result.export('anim1').set('title', 'on');
model.result.export('anim1').set('legend', 'on');
model.result.export('anim1').set('logo', 'off');
model.result.export('anim1').set('options', 'off');
model.result.export('anim1').set('fontsize', '9');
model.result.export('anim1').set('customcolor', [1 1 1]);
model.result.export('anim1').set('background', 'color');
model.result.export('anim1').set('axisorientation', 'on');
model.result.export('anim1').set('grid', 'on');
model.result.export('anim1').set('axes', 'on');
model.result.export('anim1').set('showgrid', 'on');

model.label('ISBTM.mph');

out = model;
