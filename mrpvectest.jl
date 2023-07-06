include("test.jl")

#i = Int(ceil(rand()*length(allu)))
i = 1
theta = 0 #allphi[i]
#theta = allphi[i]*rand()
mrptest = allcirclecentre[i] + sin(theta)*allu[i] + cos(theta)*allv[i]

xitest = (-1)^(i-1)*allloc_normed[sourceindices[cld(i,2)]]
xtest = allloc_normed[destindices[cld(i,2)]]


mrp_rot = MRP(mrptest)

rotatedxi = mrp_rot*xitest
xtest

isapprox(xtest,rotatedxi)

xiorig = (-1)^(i-1)*allloc[sourceindices[cld(i,2)]]
xorig = allloc[destindices[cld(i,2)]]
testrdiff = rdiff[cld(i,2)]

rotatedxiorig = mrp_rot*xiorig
rotateddist = norm(rotatedxiorig-xorig)
isapprox(rotateddist^2,testrdiff^2)

norm(xorig)
epsilon = mrptol
radialmrpvec = normalize(sin(theta)*allu[i] + cos(theta)*allv[i])
normalmrpvec = allcirclenormal[i]
phicurr = allphi[i]

randangle = pi*(2rand()-1)
mrp_rot_mod = MRP(mrptest+epsilon*radialmrpvec)
mrp_rot_mod2 = MRP(mrptest+epsilon*normalmrpvec)
mrp_rot_generror = MRP(mrptest + epsilon*cos(randangle)*radialmrpvec + epsilon*sin(randangle)*normalmrpvec)

rotatedxiorig_mod = mrp_rot_mod*xiorig
rotatedxiorig_mod2 = mrp_rot_mod2*xiorig
rotatedxiorig_generror = mrp_rot_generror*xiorig
sqdist_mod = norm(rotatedxiorig_mod-xorig)^2
sqdist_mod2 = norm(rotatedxiorig_mod2-xorig)^2
sqdist_generror = norm(rotatedxiorig_generror-xorig)^2
theoretical_sqdiff = 4*epsilon^2*sin(phicurr)^4/(1-cos(theta)*cos(phicurr))^2*norm(xorig)^2 + testrdiff^2
theoretical_max_sqdiff = 16*epsilon^2*norm(xorig)^2 + testrdiff^2
theoretical_min_sqdiff = 4*epsilon^2*norm(xorig)^2 + testrdiff^2

sqrt(theoretical_min_sqdiff)
sqrt(theoretical_sqdiff)

theoretical_sqdiff - sqdist_mod
theoretical_sqdiff - sqdist_mod2
theoretical_sqdiff - sqdist_generror

figure()
thetas_test = collect((-1:0.01:1)*phicurr)
plot(thetas_test,(4*sin(phicurr)^4)./(1.0.-cos.(thetas_test)*cos(phicurr)).^2)
display(gcf())

