import json
data = {
    'completed_in': 0.074,
    'max_id': 264043230692245504,
    'max_id_str': '264043230692245504',
    'next_page': '?page=2&max_id=264043230692245504&q=python&rpp=5',
    'page': 1,
    'query': 'python',
    'refresh_url': '?since_id=264043230692245504&q=python',
    'results': [
        {'created_at': 'Thu, 01 Nov 2012 16:36:26 +0000',
         'from_user': 'cat',
        },
        {'created_at': 'Thu, 01 Nov 2012 16:36:14 +0000',
         'from_user': 'cat',
        },
        {'created_at': 'Thu, 01 Nov 2012 16:36:13 +0000',
         'from_user': 'cat',
        },
        {'created_at': 'Thu, 01 Nov 2012 16:36:07 +0000',
         'from_user': 'cat',
        },
        {'created_at': 'Thu, 01 Nov 2012 16:36:04 +0000',
         'from_user': 'cat',
        },
    ],
    'results_per_page': 5,
    'since_id': 0,
    'since_id_str': '0'
}

#print(json.dumps(data, indent=4))
with open('f.json', 'wb') as f:
    json.dump(data, f, indent=4)
#parsed_json = json.loads(json_string)

with open('f.json', 'r') as f2:
    data2 = json.load(f2)
print(data2)
