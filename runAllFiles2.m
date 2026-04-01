
for sub = 2:8
    fprintf('%s. %d. %s. %s ...\n','prfSampleModel',sub,'ma','contour');
    prfSampleModel_symmetry(sub,8, 'ma', 'contour')
end


for sub = 1:8
    fprintf('%s. %d. %s. %s ...\n','prfSampleModel',sub,'ma','medialAxis');
    prfSampleModel_symmetry(sub,8, 'ma', 'medialAxis')
end


sym_type = {'mir', 'tap', 'par'};
for type = 1:3

    for sub = 1:8
        fprintf('%s. %d. %s. %s ...\n','prfSampleModel',sub,sym_type{type},'contour');
        prfSampleModel_symmetry(sub,8, sym_type{type}, 'contour')
        fprintf('%s. %d. %s. %s ...\n','prfSampleModel',sub,sym_type{type},'medialAxis');
        prfSampleModel_symmetry(sub,8, sym_type{type}, 'medialAxis')
        fprintf('%s. %d. %s. %s ...\n','prfSampleModel',sub,sym_type{type},'area');
        prfSampleModel_symmetry(sub,8, sym_type{type}, 'area')
    end
end