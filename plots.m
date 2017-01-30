function plots(h,kj,Shatm1,Shatm2,Shatm3,Shatm4,Shatm5,Shatm6,Shatm7,Shatm8,Shatm9)
figure;
subplot(5,2,1);
plot(h,kj,'r');
title(''),xlabel('iterations'),ylabel('kappa');
%hold on;
subplot(5,2,2);
plot(h,Shatm1,'b');
title(''),xlabel('iterations'),ylabel('Shat1');
subplot(5,2,3);
plot(h,Shatm2,'g');
title(''),xlabel('iterations'),ylabel('Shat2');
subplot(5,2,4);
plot(h,Shatm3,'c');
title(''),xlabel('iterations'),ylabel('Shat3');
subplot(5,2,5);
plot(h,Shatm4,'k');
title(''),xlabel('iterations'),ylabel('Shat4');
subplot(5,2,6);
plot(h,Shatm5,'m');
title(''),xlabel('iterations'),ylabel('Shat5');
subplot(5,2,7);
plot(h,Shatm6);
title(''),xlabel('iterations'),ylabel('Shat6');
subplot(5,2,8);
plot(h,Shatm7);
title(''),xlabel('iterations'),ylabel('Shat7');
subplot(5,2,9);
plot(h,Shatm8);
title(''),xlabel('iterations'),ylabel('Shat8');
subplot(5,2,10);
plot(h,Shatm9);
%hold off; 
title(''),xlabel('iterations'),ylabel('Shat9');

end