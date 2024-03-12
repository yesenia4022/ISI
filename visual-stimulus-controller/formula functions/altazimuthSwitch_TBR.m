function altazimuth2 = altazimuthSwitch_TBR(ori2)

if ori2 == 0 | ori2 == 180
    altazimuth2 = 'azimuth';
elseif ori2 == 90 | ori2 == 270
    altazimuth2 = 'altitude';
end 
