class NamesStorage:
    def __init__(self):
        self.loaded_names = {}

    #def __contains__(self, name):
        #"""finds out if the approximate key is in the loaded_names dictionary"""
        ##name = (vector_size, subcase_id, result_type, label, min_value, max_value)
        #key = name[:4]
        #value = name[4:]
        #assert len(key) == 4, key
        #assert len(value) == 2, value
        #if key in self.loaded_names:
            #value2 = self.loaded_names[key]
            #if value == value2:
                #return True
        #return False
        #return  in self.loaded_names

    def add(self, name):
        """
        adds the approximate name and value to the loaded_names dictionary
        """
        #print('name =', name)
        key = name[:4]
        value = name[4:]
        assert len(key) == 4, key
        #assert len(value) == 2, value
        assert key not in self.loaded_names
        self.loaded_names[key] = value
        #self.loaded_names.add(name)

    def remove(self, name):
        """removes the approximate name from the loaded_names dictionary"""
        key = name[:4]
        del self.loaded_names[key]

    def get_name_string(self, name):
        """Gets the approximate name as a string"""
        key = name[:4]
        value = name[4:]
        return '_'.join([str(k) for k in key])

    def has_close_name(self, name):
        """checks to see if the approximate key is in loaded_names"""
        key = name[:4]
        return key in self.loaded_names

    def has_exact_name(self, name):
        """
        checks to see if the approximate key is in loaded_names
        and the value is the expected value
        """
        key = name[:4]
        value = name[4:]
        return key in self.loaded_names and self.loaded_names[key] == value

