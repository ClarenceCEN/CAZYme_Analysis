

def find_value_recursive(cazy_dict,cazyme):
    global flag
    global levels
    flag = 0
    levels = cazyme
    for key,value in cazy_dict.items():
        if isinstance(value,dict):
            global n
            n += 1
            find_value_recursive(value, cazyme)
            if flag == 1:
                n -= 1
                levels =  'L'+str(n)+'_'+key+';'+levels
                return levels
        elif cazyme == value:
            n += 1
            flag = 1
            levels =  'L'+str(n)+'_'+key+';'+str(n+1)+levels
            return





def main():

    CAZYme_family = {'CAZyme':{'Enzyme_Classes':{'GlycosylTransferases':'GT',
              'GlycosideHydrolases':'GH',
              'CarbohydrateEsterases':'CE',
              'PolysaccharideLyases':'PL',
              'AuxiliaryActivities':'AA'},
              'Associated_Modules':{'Carbohydrate-BindingModules':'CBM'}}}

    global n
    n = 0

    a  = find_value_recursive(CAZYme_family,'GH')
    print(a)

if __name__  == '__main__':
    main()