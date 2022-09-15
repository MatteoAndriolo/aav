import imp
from urllib.request import urlopen
import requests
import re
from bs4 import BeautifulSoup as bs


query_url='''https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22seqmotif%22%2C%22parameters%22%3A%7B%22value%22%3A%22AWDAWDAWD%22%2C%22pattern_type%22%3A%22simple%22%2C%22sequence_type%22%3A%22protein%22%7D%7D%5D%7D%2C%22return_type%22%3A%22entry%22%2C%22request_options%22%3A%7B%22paginate%22%3A%7B%22start%22%3A0%2C%22rows%22%3A25%7D%2C%22results_content_type%22%3A%5B%22experimental%22%5D%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22desc%22%7D%5D%2C%22scoring_strategy%22%3A%22combined%22%7D%2C%22request_info%22%3A%7B%22query_id%22%3A%22169051ad6768f923fe9422feea59d787%22%7D%7D'''
search_motif="MATTALA"

url="https://www.rcsb.org/search?request="
query:dict={
  "query": {
    "type": "terminal",
    "service": "seqmotif",
    "parameters": {
      "value": search_motif,
      "pattern_type": "prosite",
      "sequence_type": "protein"
    }
  },
  "return_type": "polymer_entity"
}
url=url+str(query)
url=url.replace(" ","").replace("'",'"')

# selenium
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support import expected_conditions as EC

def _generate_url_searchmotif(motif:str)->str:
    url="https://www.rcsb.org/search?request="
    query:dict={
            "query": {
                "type": "terminal",
                "service": "seqmotif",
                "parameters": {
                "value": motif,
                "pattern_type": "prosite",
                "sequence_type": "protein"
                }
            },
        "return_type": "polymer_entity"
        }
    url=url+str(query)
    url=url.replace(" ","").replace("'",'"')
    return url

def _generate_url_structure(pdb_id:str):
  url=f"https://www.rcsb.org/structure/{pdb_id}"
  return url

options = Options()
options.binary_location = "/usr/sbin/google-chrome-stable"
options.add_argument("user-data-dir=/home/matteo/cache")
driver = webdriver.Chrome(chrome_options = options, executable_path='/home/matteo/chromedriver_linux64/chromedriver')
driver.get(url)
print("Chrome Browser Invoked")

#WebDriverWait wait = new WebDriverWait(webDriver, timeoutInSeconds);
#wait.until(ExpectedConditions.visibilityOfElementLocated(By.id<locator>));
try:
  element = WebDriverWait(driver, 20).until(
    EC.presence_of_element_located((By.CLASS_NAME, "results-item"))
  )
except NoSuchElementException:
  print("not found")

##driver.find_element((By.CLASS_NAME, "results-item"))
for el in driver.find_elements(By.TAG_NAME,"h3"):
  print(re.findall(r"[A-Z\d]{4}",el.accessible_name)[0])
#print(driver.page_source.count("results-item"))
try:
    WebDriverWait(driver, 20).until(
        EC.presence_of_element_located((By.CLASS_NAME, "results-item"))
    )
    elements=driver.find_elements(By.TAG_NAME,"h3")
    pdb_id=re.findall(r"[A-Z\d]{4}",elements[0].accessible_name)[0]
except NoSuchElementException:
    print("not found")
driver.get(_generate_url_structure(pdb_id))
input("wait press button")

driver.quit()



## Using requests
#print(url.replace(" ","").replace("'",'"'))
#print(50*"#")
#response=requests.get(url+str(query))
#print(response.status_code)
#soup = bs(response.text, 'html.parser')
#print(soup.prettify())
#print(soup.find_all("div",class_="col-md-9 col-xs-12 results-item-info"))

#page=urlopen(url)
#html=page.read().decode("utf-8")
#soup = bs(html, 'html.parser')
#print(soup.prettify())
#print(soup.find_all("div",class_="col-md-9 col-xs-12 results-item-info"))