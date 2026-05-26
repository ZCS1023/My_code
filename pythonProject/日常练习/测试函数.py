import unittest
'''from 拼接名字 import get_name

class NameTest(unittest.TestCase):
    def test1(self):
        full_name=get_name('James','Lebron')
        self.assertEqual(full_name,'James Lebron')
    def test2(self):
        full_name=get_name('Feng','Zhang','Li')
        self.assertEqual(full_name,'Feng Li Zhang')

if __name__=='__main__':
    unittest.main()'''

from 城市与国家 import city_country
class test_city_country(unittest.TestCase):
    def test(self):
        full_name=city_country('Beijing','China')
        self.assertEqual(full_name,'Beijing,China')

if __name__=='__main__':
    unittest.main()