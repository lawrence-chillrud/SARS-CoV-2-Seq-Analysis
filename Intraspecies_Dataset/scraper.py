import selenium as se
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
import traceback
import logging
import os
import time

# 274 in total (b/w the two lists)... 
UNFOUND_1 = ['MT326187', 'MT326190', 'NC_045512', 'MT339041', 'MT334529', 'MT334538', 'MT334542', 'MT334544', 'MT334547', 'MT326178', 'MT324062', 'MT328032', 'MT328035', 'MT334533', 'MT334541', 'MT334523', 'MT334537', 'MT334540', 'MT334535', 'MT334536', 'MT334546', 'MT326167', 'MT339039', 'MT334526', 'MT293204', 'MT334534', 'MT326129', 'MT326052', 'MT326065', 'MT326117', 'MT326092', 'MT293159', 'MT324684', 'MT326132', 'MT293196', 'MT293213', 'MT326106', 'MT293219', 'MT328034', 'MT325563', 'MT325565', 'MT325566', 'MT325568', 'MT325569', 'MT325590', 'MT325640', 'MT325606', 'MT325607', 'MT325608', 'MT325609', 'MT325610', 'MT325611', 'MT325616', 'MT325618', 'MT325619', 'MT325620', 'MT325622', 'MT325623', 'MT325624', 'MT325599', 'MT325600', 'MT325601', 'MT325602', 'MT325612', 'MT325613', 'MT325615', 'MT325617', 'MT325625', 'MT325573', 'MT325574', 'MT325577', 'MT325579', 'MT325586', 'MT325592', 'MT325593', 'MT325594', 'MT325598', 'MT325605', 'MT325626', 'MT325627', 'MT325633', 'MT325634', 'MT326069', 'MT325561', 'MT325571', 'MT325572', 'MT325575', 'MT325583', 'MT325587', 'MT325588', 'MT325589', 'MT325596', 'MT325597', 'MT325603', 'MT325604', 'MT325614', 'MT325621', 'MT325629', 'MT325630', 'MT325631', 'MT325632', 'MT325635', 'MT325636', 'MT325637', 'MT325638', 'MT325639', 'MT325564', 'MT325567', 'MT325584', 'MT325585', 'MT325562', 'MT325570', 'MT325576', 'MT325578', 'MT325580', 'MT325581', 'MT325582', 'MT325591', 'MT325595', 'MT325628', 'MN997409', 'MT293198', 'MN988668', 'MN988669', 'MT326095', 'MT334561', 'MT326150', 'MT326096', 'MT328033', 'MT326153', 'MT326100', 'MT326063'] 

UNFOUND_2 = ['MT326088', 'MT326118', 'MT326048', 'MT326066', 'MT326067', 'MT326164', 'MT293218', 'MT066175', 'MT326061', 'MT293220', 'MT293187', 'MT246469', 'MT326140', 'MT326103', 'MT326042', 'MT322398', 'MT263432', 'MT326116', 'MT326056', 'MT293186', 'MT293225', 'MT334531', 'MT326134', 'MT326068', 'MT326027', 'MT326154', 'MT326029', 'MT334539', 'MT326137', 'MT326127', 'MT326162', 'MT326099', 'MT326055', 'MT326097', 'MT293216', 'MT326023', 'MT326113', 'MT326086', 'MT326119', 'MT012098', 'MT326051', 'MT326087', 'MT326030', 'MT293200', 'MT326089', 'MT326120', 'MT293171', 'MT293162', 'MT293224', 'MT326071', 'MT246482', 'MT326105', 'MT293222', 'MT326180', 'MT263400', 'MT263425', 'MT334557', 'MT326076', 'MT326074', 'MT293165', 'MT326078', 'MT326128', 'MT326159', 'MT326075', 'MT326081', 'MT188341', 'MT326070', 'MT326093', 'MT326082', 'MT326049', 'MT326191', 'MT327745', 'MT326133', 'MT259285', 'MT326169', 'MT326058', 'MT293190', 'MT263412', 'MT326125', 'MT326084', 'MT326028', 'MT326185', 'MT291833', 'MT326040', 'MT326148', 'MT326112', 'MT326104', 'MT263386', 'MT339040', 'MT291826', 'MT293167', 'MT326085', 'MT326135', 'MT326090', 'MT263441', 'MT326039', 'MT326130', 'MT293169', 'MT263434', 'MT292572', 'MT292574', 'MT326035', 'MT293188', 'MT326111', 'MT326046', 'MT293175', 'MT259266', 'MT263456', 'MT246485', 'MT326041', 'MT293164', 'MT326138', 'MT326124', 'MT326160', 'MT326053', 'MT326166', 'MT326152', 'MT293184', 'MT326094', 'MT293158', 'MT326036', 'MT326151', 'MT326123', 'MT326158', 'MT326184', 'MT326149', 'MT326054', 'MT293191', 'MT326091', 'MT326136', 'MT326024', 'MT326083', 'MT293199', 'MT326121', 'MT326181', 'MT293185', 'MT293197', 'MT293163', 'MT326161', 'MT326173', 'MT293180', 'MT259284']

def extractSeqs():
    accessions = []
    f = open("sequences.fasta.txt", "r")
    for line in f:
        if line[0] == ">":
            if line[1] == "N":
                accessions.append(line[1:10])
            else:
                accessions.append(line[1:9])

    return accessions

def search(accs):
    os.chdir('../../../../../../Downloads')
    print(os.getcwd())
    found, unfound = [], []
    total_found = 0
    total_unfound = 0
    for acc in accs:
        driver = se.webdriver.Chrome('/Users/lawrencechillrud/Desktop/chromedriver')
        driver.get('https://www.ncbi.nlm.nih.gov/')
        search_box = driver.find_element_by_name('term')
        search_box.send_keys(acc)
        search_button = driver.find_element_by_id('search')
        search_button.click()
        # find and click on link:
        try:
            link = driver.find_element_by_id('sequence_title')
            link.click()
            maincontent = driver.find_element_by_id('maincontent')
            content = maincontent.find_element_by_class_name('content')
            sub_content = content.find_element_by_id('viewercontent1')
            driver.implicitly_wait(3)
            seq = sub_content.find_element_by_class_name('sequence')
            div = content.find_element_by_class_name('rprtheader')
            acc_complete = div.find_element_by_class_name('itemid').text[9:]

            # uncomment if doesn't work and delete all of 'testing stuff here' below:
            #cds_id = "feature_%s_CDS_1" % (acc_complete) # CHANGE 'gene' BACK TO 'CDS' if broken
            #span = seq.find_element_by_id(cds_id)
            # testing stuff here
            n = 1
            for i in range(4):
                cds_id = "feature_%s_CDS_%d" % (acc_complete, i) # CHANGE 'gene' BACK TO 'CDS' if broken
                span = seq.find_element_by_id(cds_id)
                #print(span.text)
                if ('gene="S"' in span.text) or ("S protein" in span.text) or ("surface glycoprotein" in span.text) or ("spike" in span.text):
                    n = i
                    break

            # end testing
            #print(span.text)
            final_link = span.find_element_by_css_selector('a').get_attribute('href')
            driver.get(final_link)

            # now just need to download the file
            mc = driver.find_element_by_id('maincontent')
            c = mc.find_element_by_class_name('content')
            d = c.find_element_by_id('seqsendto')
            driver.execute_script("arguments[0].setAttribute('aria-expanded','true')", d)
            d2 = c.find_element_by_id('send_to_topmenu')
            driver.execute_script("arguments[0].setAttribute('aria-hidden','false')", d2)
            driver.execute_script("arguments[0].setAttribute('style','top: 98.2812px; left: 211.656px; display: block;')", d2)
            d3 = c.find_element_by_id('submenu_complete_rec')
            driver.execute_script("arguments[0].setAttribute('style','')", d3)

            # now click radio button:
            radio = d3.find_element_by_id('dest_File')
            radio.click()
            last_section = d3.find_element_by_id('file_format')
            for option in last_section.find_elements_by_tag_name('option'):
                if option.text == 'GenBank':
                    option.click()
                    break

            submenu_file = d2.find_element_by_id('submenu_File')
            lastb = submenu_file.find_element_by_tag_name('button')
            time.sleep(1)
            lastb.click()
            time.sleep(8)
            # move to correct place:
            new_file_name = '../Documents/CU_SeniorYear/Spring2020/ComputationalGenomics/Final_Project/SARS-CoV-2-seqs/scraped/%s.txt' % (acc)
            os.rename('sequence.gb', new_file_name)
            print("Completed: %s" % (acc_complete))
            found.append(acc)
            total_found += 1
            driver.quit()
        except Exception as e:
            print("Error processing: %s" % (acc))
            logging.error(traceback.format_exc())
            # Logs the error appropriately. 
            unfound.append(acc)
            total_unfound += 1
            driver.quit()

    return found, unfound, total_found, total_unfound

if __name__ == "__main__":
    #accessionsIwant = extractSeqs()
    #print(accessionsIwant)
    f, u, tf, tu = search(UNFOUND_1)
    print("Total found: %d / %d" % (tf, tf+tu))
    print(f)
    print("Total unfound: %d / %d" % (tu, tu+tf))
    print(u)
    #print(accessionsIwant[-100:])
