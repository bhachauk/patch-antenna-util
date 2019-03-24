function toggle(evt, id) {
    makeDisplay ('show', 'none')
    x = document.getElementsByClassName('currentBtn');
    for (i = 0; i < x.length; i++) {
        x[i].className = x[i].className.replace("currentBtn", "cmnBtn");
    }
    evt.currentTarget.className = "currentBtn";
    makeDisplayId(id, 'block');
}

function makeDisplay (className, disType)
{
    var list = document.getElementsByClassName(className);
     for (i = 0; i < list.length; i++) {
       list[i].style.display = disType;
     }
}

function makeDisplayId (id, disType)
{
    document.getElementById (id).style.display = disType;
}

function simulate()
{
    makeDisplayId('totalContainer','block');
}